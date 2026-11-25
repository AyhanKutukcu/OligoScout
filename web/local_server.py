"""Local OligoScout HTTP application. Run from WSL with Python 3.10+."""
from __future__ import annotations

import argparse
import concurrent.futures
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import os
from pathlib import Path
import re
import tempfile
import threading
import traceback
from urllib.parse import urlsplit
import uuid

from local_engine import LocalEngine, clean_json

WEB = Path(__file__).resolve().parent
JOB_ID = re.compile(r"[a-f0-9]{24}")
# Page addresses of the single-page app (oligoscout.org/design, /docs/pair/...):
# each is answered with index.html and the browser-side router picks the view.
APP_ROUTE = re.compile(r"/(?:overview|design|multiplex|pair|analysis|catalog|docs|contact|"
                       r"panel|probes|melt|validation|"
                       # older Turkish addresses, rewritten to the English ones in the browser
                       r"genel|tasarla|multipleks|cift|yerel|katalog|dokumantasyon|iletisim|"
                       r"prob|erime|dogrulama)(?:/[A-Za-z0-9_.:%-]+){0,3}/?")
MAX_BODY = 64000
MAX_JOBS = 40
MAX_PENDING = 4
PAIR_TOOLS = {"qc", "variants", "melt", "probe", "isoform"}
TOOLS = PAIR_TOOLS | {"design", "multiplex_design", "scan", "thermo", "blocker", "vpcr"}

# --- Public deployment (reverse proxy) ---------------------------------------
# Comma-separated hostnames allowed behind a TLS-terminating reverse proxy,
# e.g. OLIGOSCOUT_ALLOWED_HOSTS="oligoscout.org,www.oligoscout.org".
# Empty (default) keeps the original localhost-only behaviour unchanged.
PUBLIC_HOSTS = tuple(h.strip().lower() for h in os.getenv("OLIGOSCOUT_ALLOWED_HOSTS", "").split(",") if h.strip())
PUBLIC_SCHEME = os.getenv("OLIGOSCOUT_PUBLIC_SCHEME", "https").strip().lower()


class RequestError(ValueError):
    def __init__(self, message, status=HTTPStatus.BAD_REQUEST):
        super().__init__(message)
        self.status = status


class JobQueue:
    def __init__(self, motor):
        self.motor = motor
        self.jobs = {}
        self.lock = threading.Lock()
        self.worker = concurrent.futures.ThreadPoolExecutor(max_workers=1, thread_name_prefix="oligoscout")

    def get(self, key):
        with self.lock:
            job = self.jobs.get(key)
            if job is None:
                raise RequestError("Analiz bu oturumda bulunamadı.", HTTPStatus.NOT_FOUND)
            return dict(job)

    def submit(self, data):
        tool = data.get("tool")
        if not isinstance(tool, str) or tool not in TOOLS:
            raise RequestError("Geçerli bir analiz seçin.")
        source = None
        selected = None
        if tool in PAIR_TOOLS or tool == "scan":
            source_id = data.get("source", "project")
            if source_id == "project":
                source = self.motor.project()
            elif isinstance(source_id, str) and JOB_ID.fullmatch(source_id):
                design_job = self.get(source_id)
                if design_job["tool"] not in ("design", "multiplex_design") or design_job["status"] != "done":
                    raise RequestError("Önce primer tasarımının tamamlanmasını bekleyin.")
                source = design_job["result"]
            else:
                raise RequestError("Analiz kaynağı geçersiz.")
            lookup = {p["id"]: p for p in source["pairs"]}
            ids = data.get("pairs")
            if (not isinstance(ids, list) or not 1 <= len(ids) <= 50
                    or any(not isinstance(key, str) or key not in lookup for key in ids)
                    or len(set(ids)) != len(ids)):
                raise RequestError("Panelden geçerli primer çifti veya çiftleri seçin.")
            if tool != "scan" and len(ids) != 1:
                raise RequestError("Bu analiz için bir primer çifti seçin.")
            selected = [lookup[key] for key in ids]
        job_id = uuid.uuid4().hex[:24]
        with self.lock:
            if sum(x["status"] in ("queued", "running") for x in self.jobs.values()) >= MAX_PENDING:
                raise RequestError("Analiz sırası dolu. Çalışan bir analiz tamamlanınca tekrar deneyin.",
                                   HTTPStatus.TOO_MANY_REQUESTS)
            if len(self.jobs) >= MAX_JOBS:
                oldest = next(key for key, job in self.jobs.items() if job["status"] in ("done", "error"))
                del self.jobs[oldest]
            self.jobs[job_id] = {"id": job_id, "tool": tool, "status": "queued", "stage": "Sırada"}
        self.worker.submit(self.run, job_id, tool, data, source, selected)
        return self.get(job_id)

    def update(self, key, **fields):
        with self.lock:
            self.jobs[key].update(fields)

    def run(self, key, tool, data, source, selected):
        def stage(name):
            self.update(key, status="running", stage=name)

        stage("Yerel motor hazırlanıyor")
        try:
            with tempfile.TemporaryDirectory(prefix="oligoscout-web-") as tmp:
                work = Path(tmp)
                if tool == "design":
                    result = self.motor.design(data, stage, work)
                elif tool == "multiplex_design":
                    from multiplex_design import run
                    result = run(self.motor, data, stage, work)
                elif tool == "scan":
                    result = self.motor.scan(source, selected, work, stage, data)
                elif tool in PAIR_TOOLS:
                    result = self.motor.single(tool, source, selected[0], data)
                else:
                    result = self.motor.oligo(tool, data, work)
            self.update(key, status="done", stage="Tamamlandı", result=clean_json(result))
        except ValueError as exc:
            self.update(key, status="error", stage="Analiz tamamlanamadı", error=str(exc))
        except Exception:
            traceback.print_exc()
            self.update(key, status="error", stage="Analiz tamamlanamadı",
                        error="Yerel motor hatası. Ayrıntı için sunucu terminaline bakın.")


class Handler(BaseHTTPRequestHandler):
    server_version = "OligoScout"

    def setup(self):
        super().setup()
        self.connection.settimeout(15)

    def log_message(self, format, *args):
        if args and isinstance(args[0], str) and args[0].startswith(("GET / ", "POST ")):
            super().log_message(format, *args)

    def send_data(self, status, payload, content_type="application/json; charset=utf-8"):
        body = payload if isinstance(payload, bytes) else json.dumps(
            payload, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode()
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Content-Security-Policy",
                         "default-src 'none'; script-src 'unsafe-inline'; style-src 'unsafe-inline'; "
                         "img-src 'self' data:; connect-src 'self'; font-src 'self' data:; "
                         "base-uri 'none'; frame-ancestors 'none'; form-action 'self'; frame-src https://www.google.com/maps/embed")
        self.end_headers()
        self.wfile.write(body)

    def check_host(self):
        host = self.headers.get("Host", "")
        port = self.server.server_port
        local = (f"127.0.0.1:{port}", f"localhost:{port}")
        if host.lower() in PUBLIC_HOSTS:
            # Served publicly behind a reverse proxy that preserves the Host header.
            origin = self.headers.get("Origin")
            if origin is not None and origin != f"{PUBLIC_SCHEME}://{host.lower()}":
                raise RequestError("İstek bu uygulamadan gelmeli.", HTTPStatus.FORBIDDEN)
        elif host in local:
            origin = self.headers.get("Origin")
            if origin is not None and origin != "http://" + host:
                raise RequestError("İstek bu yerel uygulamadan gelmeli.", HTTPStatus.FORBIDDEN)
        else:
            raise RequestError("Yerel sunucuyu localhost üzerinden açın.", HTTPStatus.FORBIDDEN)
        if self.headers.get("Sec-Fetch-Site") == "cross-site" and self.command == "POST":
            raise RequestError("Çapraz site isteği kabul edilmiyor.", HTTPStatus.FORBIDDEN)

    def handle_error(self, exc):
        code = exc.status if isinstance(exc, RequestError) else HTTPStatus.BAD_REQUEST
        self.send_data(code, {"error": str(exc)})

    def do_GET(self):
        try:
            self.check_host()
            path = urlsplit(self.path).path
            if path in ("/", "/index.html") or APP_ROUTE.fullmatch(path):
                self.send_data(HTTPStatus.OK, (WEB / "index.html").read_bytes(), "text/html; charset=utf-8")
            elif path == "/api/status":
                self.send_data(HTTPStatus.OK, self.server.jobs.motor.status())
            elif path == "/api/project":
                self.send_data(HTTPStatus.OK, self.server.jobs.motor.project())
            elif path.startswith("/api/jobs/") and JOB_ID.fullmatch(path.rsplit("/", 1)[-1]):
                self.send_data(HTTPStatus.OK, self.server.jobs.get(path.rsplit("/", 1)[-1]))
            else:
                self.send_data(HTTPStatus.NOT_FOUND, {"error": "Sayfa bulunamadı."})
        except (RequestError, ValueError) as exc:
            self.handle_error(exc)

    def do_POST(self):
        try:
            self.check_host()
            if urlsplit(self.path).path != "/api/jobs":
                raise RequestError("İşlem bulunamadı.", HTTPStatus.NOT_FOUND)
            if self.headers.get_content_type() != "application/json":
                raise RequestError("İstek JSON biçiminde olmalı.", HTTPStatus.UNSUPPORTED_MEDIA_TYPE)
            size = int(self.headers.get("Content-Length", "0"))
            if not 1 <= size <= MAX_BODY:
                raise RequestError("İstek boş veya çok büyük.", HTTPStatus.REQUEST_ENTITY_TOO_LARGE)
            body = json.loads(self.rfile.read(size))
            if not isinstance(body, dict):
                raise RequestError("İstek bir JSON nesnesi olmalı.")
            self.send_data(HTTPStatus.ACCEPTED, self.server.jobs.submit(body))
        except (RequestError, ValueError) as exc:
            self.handle_error(exc)


def main():
    parser = argparse.ArgumentParser(description="OligoScout yerel web sunucusu")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--host", default="127.0.0.1",
                        help="0.0.0.0 only inside a container or behind a reverse proxy")
    args = parser.parse_args()
    if not (WEB / "index.html").is_file():
        parser.error("Önce python3 web/build_web_data.py ve python3 web/build_site.py çalıştırın.")
    with ThreadingHTTPServer((args.host, args.port), Handler) as server:
        server.jobs = JobQueue(LocalEngine())
        print(f"OligoScout hazır. http://127.0.0.1:{args.port}/#genel", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            print("\nYerel sunucu durduruluyor.", flush=True)
        finally:
            server.jobs.worker.shutdown(wait=True, cancel_futures=True)


if __name__ == "__main__":
    main()
