// OligoScout tarayıcı motoru — reference_py/oligoscout_thermo'nun birebir portu.
// Parametre tabloları kopyalanmaz: veri paketindeki engine_params (params.py dışa aktarımı) kullanılır.
const OS = (() => {
  let WC=null, INIT=null, SYM=[0,-1.4], MM=null, LOOP=null, ATB=0.5, R=1.9872, K0=273.15;
  const T37=310.15;
  const COMP={A:'T',T:'A',G:'C',C:'G',N:'N'};
  const PAIR=(x,y)=>(x==='A'&&y==='T')||(x==='T'&&y==='A')||(x==='G'&&y==='C')||(x==='C'&&y==='G');
  const rc=s=>{let o='';for(let i=s.length-1;i>=0;i--)o+=COMP[s[i]]||'N';return o;};
  const gcf=s=>{let n=0;for(const c of s)if(c==='G'||c==='C')n++;return s.length?n/s.length:0;};
  function setParams(p){WC=p.wc.nn;INIT=p.wc.init;SYM=p.wc.symmetry||SYM;MM=p.mismatch;LOOP=p.loop;ATB=p.at_bulge_closure_dg37;R=p.R_cal||R;K0=p.kelvin||K0;}
  const step=s=>WC[s]||WC[rc(s)];
  function stackE(seg){let h=0,s=0;for(let i=0;i<seg.length-1;i++){const v=step(seg.substr(i,2));h+=v[0];s+=v[1];}return [h,s];}
  const initE=b=>INIT[(b==='A'||b==='T')?'AT':'GC'];
  function helixE(seg){let [h,s]=stackE(seg);for(const b of [seg[0],seg[seg.length-1]]){const v=initE(b);h+=v[0];s+=v[1];}return [h,s];}
  const selfc=s=>s===rc(s);
  function watsonCrick(seq){let [h,s]=helixE(seq);const sc=selfc(seq);if(sc){h+=SYM[0];s+=SYM[1];}return {h,s,selfc:sc};}
  // ---- tuz (salt.py) ----
  const _A=3.92e-5,_B=-9.11e-6,_C=6.26e-5,_D=1.42e-5,_E=-4.82e-4,_F=5.25e-4,_G=8.31e-5;
  function ow04K(tmK,na,f){const l=Math.log(na);return 1/(1/tmK+(4.29*f-3.95)*1e-5*l+9.40e-6*l*l);}
  function ow08K(tmK,mg,mon,f,n){
    if(mg<=0) return [ow04K(tmK,mon,f),'monovalent_only'];
    const ratio=mon>0?Math.sqrt(mg)/mon:Infinity;
    if(ratio<0.22) return [ow04K(tmK,mon,f),'monovalent_dominant'];
    let a=_A,d=_D,g=_G,reg='divalent_dominant';
    if(ratio<6){const l=Math.log(mon);a=_A*(0.843-0.352*Math.sqrt(mon)*l);d=_D*(1.279-4.03e-3*l-8.03e-3*l*l);g=_G*(0.486-0.258*l+5.25e-3*l*l*l);reg='mixed_mono_divalent';}
    const lm=Math.log(mg);
    return [1/(1/tmK+a+_B*lm+f*(_C+d*lm)+(1/(2*(n-1)))*(_E+_F*lm+g*lm*lm)),reg];
  }
  const naEqVA=(mon,mg,dntp)=>mon+120*Math.sqrt(Math.max(0,(mg-dntp)*1e3))/1e3;
  // c: {na,mg,dntp (mM), pnm (nM), ta (°C)}
  function M(c){return {mon:c.na/1e3,mg:c.mg/1e3,dntp:c.dntp/1e3,free:Math.max(0,(c.mg-c.dntp)/1e3),C:c.pnm*1e-9,ta:c.ta};}
  function saltTmK(tm1K,seq,c,model){const k=M(c),f=gcf(seq),n=seq.length;
    if(model==='owczarzy_2004')return [ow04K(tm1K,k.mon,f),'monovalent'];
    if(model==='von_ahsen_2001')return [ow04K(tm1K,naEqVA(k.mon,k.mg,k.dntp),f),'na_equivalent'];
    return ow08K(tm1K,k.free,k.mon,f,n);}
  // concentration.build(primer_excess)
  function thermo(seq,c,model='owczarzy_2008'){
    const k=M(c),d=watsonCrick(seq);const ceff=k.C;
    const tm1K=d.h*1000/(d.s+R*Math.log(ceff));
    const [tmK,regime]=saltTmK(tm1K,seq,c,model);
    const sEff=d.h*1000/tmK-R*Math.log(ceff);
    const dGcal=T=>d.h*1000-(T+K0)*sEff;
    const theta=T=>{const x=Math.exp(-dGcal(T)/(R*(T+K0)))*ceff;return x/(1+x);};
    return {dH:d.h,dS:d.s,dS_eff:sEff,tm1M:tm1K-K0,tm:tmK-K0,regime,theta,thetaTa:theta(k.ta),dGta:dGcal(k.ta)/1000,dG37:d.h-T37*d.s/1000,ceff};
  }
  function saltSpread(seq,c){const o={};for(const m of ['owczarzy_2004','owczarzy_2008','von_ahsen_2001'])o[m]=thermo(seq,c,m).tm;
    const used=c.mg>0?['owczarzy_2008','von_ahsen_2001']:['owczarzy_2004','owczarzy_2008'];const v=used.map(m=>o[m]);o.spread=Math.max(...v)-Math.min(...v);o.used=used;return o;}
  function dG3(seq){const [h,s]=stackE(seq.slice(-5));return h-T37*s/1000;}
  function halves(seq){const m=Math.floor(seq.length/2);const a=stackE(seq.slice(0,m)),b=stackE(seq.slice(m));return {d5:a[0]-T37*a[1]/1000,d3:b[0]-T37*b[1]/1000};}
  // ---- align.py ----
  function mmAdjust(top,bot,pos){let dh=0,ds=0;
    for(const i of [pos-1,pos]){const ts=top.substr(i,2),bs=bot.substr(i,2);
      const k1=ts+'/'+bs,k2=bs.split('').reverse().join('')+'/'+ts.split('').reverse().join('');const key=k1<k2?k1:k2;
      const w=step(ts),m=MM[key];if(!m)return null;dh+=m[0]-w[0];ds+=m[1]-w[1];}
    return [dh,ds];}
  function mkRes(o){return o;}
  function scanContig(sa,rb,mode){
    const na=sa.length,nb=rb.length;let best=null;
    for(let off=-(nb-1);off<na;off++){
      const lo=Math.max(0,off),hi=Math.min(na,off+nb);if(hi-lo<2)continue;
      const paired=[];for(let k=lo;k<hi;k++)paired.push(PAIR(sa[k],rb[k-off]));
      let start=null,prev=null;
      for(let idx=0;idx<=paired.length;idx++){
        if(idx<paired.length&&paired[idx]){if(start===null)start=idx;continue;}
        if(start===null)continue;
        const rlo=lo+start,rhi=lo+idx;start=null;const cands=[];
        if(rhi-rlo>=2){cands.push([rlo,null]);if(prev){const [plo,phi]=prev;if(rlo===phi+1&&phi-plo>=2)cands.push([plo,phi]);}}
        prev=[rlo,rhi];
        for(const [left,mm] of cands){
          const terminal=(rhi===na)||(left===off);
          if(mode==='3prime'&&!terminal)continue;
          const top=sa.slice(left,rhi);let [dH,dS]=helixE(top);
          if(mm!==null){const bot=rb.slice(left-off,rhi-off);const adj=mmAdjust(top,bot,mm-left);if(!adj)continue;dH+=adj[0];dS+=adj[1];}
          const dG=dH-310.15*dS/1000;if(best&&dG>=best.dG)continue;
          const pad=Math.max(0,-off);const tl=' '.repeat(pad)+sa,bl=' '.repeat(Math.max(0,off))+rb;const w=Math.max(tl.length,bl.length);
          const marks=Array(w).fill(' ');for(let x=pad+left;x<pad+rhi;x++)marks[x]='|';const cnt=mm!==null?1:0;if(cnt)marks[pad+mm]='X';
          best={dH,dS,dG,offset:off,helix_len:rhi-left-cnt,mode,alignment:[tl.padEnd(w),marks.join(''),bl.padEnd(w)],three_prime_paired:terminal,mismatch_count:cnt,bulge_count:0};
        }
      }
    }
    return best;
  }
  function scanBulge(sa,rb,mode){
    const na=sa.length,nb=rb.length;let best=null;const bh=0,bs=-1000*LOOP['BULGE:1'][0]/T37;
    for(let off=-(nb-1);off<na;off++){
      const lo=Math.max(0,off),hi=Math.min(na,off+nb);
      for(const strand of ['b','a']){
        const o2=strand==='b'?off-1:off+1;const lo2=Math.max(0,o2),hi2=Math.min(na,o2+nb);
        for(let h1=lo+2;h1<=hi;h1++){
          if(!PAIR(sa[h1-1],rb[h1-1-off]))continue;
          let left=h1-1;while(left-1>=lo&&PAIR(sa[left-1],rb[left-1-off]))left--;
          if(h1-left<2)continue;
          const l2=strand==='b'?h1:h1+1;if(l2<lo2||l2>=hi2)continue;
          if(!PAIR(sa[l2],rb[l2-o2]))continue;
          let h2=l2;while(h2<hi2&&PAIR(sa[h2],rb[h2-o2]))h2++;
          if(h2-l2<2)continue;
          let dH=0,dS=0;for(const seg of [sa.slice(left,h1),sa.slice(l2,h2)]){const e=stackE(seg);dH+=e[0];dS+=e[1];}
          for(const b of [sa[left],sa[h2-1]]){const e=initE(b);dH+=e[0];dS+=e[1];}
          dH+=bh;dS+=bs;
          const atc=('AT'.includes(sa[h1-1])?1:0)+('AT'.includes(sa[l2])?1:0);dS-=1000*ATB*atc/T37;
          let e;if(strand==='b')e=step(sa[h1-1]+sa[h1]);else{const j=h1-1-off;e=step(rb[j+1]+rb[j]);}
          dH+=e[0];dS+=e[1];const dG=dH-T37*dS/1000;
          const terminal=(h2===na)||(left===off);if(mode==='3prime'&&!terminal)continue;if(best&&dG>=best.dG)continue;
          const pad=Math.max(0,-off);let tl,bl,rs;
          if(strand==='b'){tl=' '.repeat(pad)+sa.slice(0,h1)+'-'+sa.slice(h1);bl=' '.repeat(Math.max(0,off))+rb;rs=1;}
          else{const m=h1-off;tl=' '.repeat(pad)+sa;bl=' '.repeat(Math.max(0,off))+rb.slice(0,m)+'-'+rb.slice(m);rs=0;}
          const w=Math.max(tl.length,bl.length);const marks=Array(w).fill(' ');
          for(let k=left;k<h1;k++)marks[pad+k]='|';for(let k=l2;k<h2;k++)marks[pad+k+rs]='|';marks[pad+h1]='B';
          best={dH,dS,dG,offset:off,helix_len:(h1-left)+(h2-l2),mode,alignment:[tl.padEnd(w),marks.join(''),bl.padEnd(w)],three_prime_paired:terminal,mismatch_count:0,bulge_count:1,bulge_strand:strand};
        }
      }
    }
    return best;
  }
  function bestDimer(a,b,mode='any',allowBulge=true){
    const rb=b.split('').reverse().join('');let best=scanContig(a,rb,mode);
    if(allowBulge){const bu=scanBulge(a,rb,mode);if(bu&&(!best||bu.dG<best.dG))best=bu;}
    return best;
  }
  // IDT tarzı çoklu yapı listesi: her kayma için en iyi kesintisiz/uyumsuz helis
  function dimerList(a,b,n=5){
    const rb=b.split('').reverse().join('');const out=[];
    for(let off=-(rb.length-1);off<a.length;off++){const r=scanContigOffset(a,rb,off);if(r)out.push(r);}
    out.sort((x,y)=>x.dG-y.dG);return out.slice(0,n);
  }
  function scanContigOffset(sa,rb,off){ // tek kayma
    const na=sa.length,nb=rb.length,lo=Math.max(0,off),hi=Math.min(na,off+nb);if(hi-lo<2)return null;
    let best=null;let start=null;
    for(let k=lo;k<=hi;k++){const p=k<hi&&PAIR(sa[k],rb[k-off]);
      if(p){if(start===null)start=k;continue;}
      if(start!==null&&k-start>=2){const top=sa.slice(start,k);const [dH,dS]=helixE(top);const dG=dH-T37*dS/1000;
        if(!best||dG<best.dG){const pad=Math.max(0,-off);const tl=' '.repeat(pad)+sa,bl=' '.repeat(Math.max(0,off))+rb;const w=Math.max(tl.length,bl.length);const m=Array(w).fill(' ');for(let x=pad+start;x<pad+k;x++)m[x]='|';
          best={dG,dH,dS,offset:off,helix_len:k-start,alignment:[tl.padEnd(w),m.join(''),bl.padEnd(w)],three_prime_paired:(k===na)||(start===off),mismatch_count:0,bulge_count:0};}}
      start=null;}
    return best;
  }
  function hairpinStem(s,minLoop=3){
    let best=null;const n=s.length;
    for(let i=0;i<n;i++)for(let j=i+minLoop+1;j<n;j++){
      let L=0;while(i-L>=0&&j+L<n&&PAIR(s[i-L],s[j+L]))L++;
      if(L<2)continue;const seg=s.slice(i-L+1,i+1);const [dH,dS]=helixE(seg);const dG=dH-310.15*dS/1000;
      if(!best||dG<best.dG){let db=Array(n).fill('.');for(let k=0;k<L;k++){db[i-k]='(';db[j+k]=')';}
        best={dG,dH,dS,i,j,helix_len:L,loop:j-i-1,db:db.join(''),three_prime_paired:(j+L-1===n-1)};}
    }
    return best;
  }
  // ---- kompozisyon (design.py / seqpair.py) ----
  const longestHomo=s=>{if(!s)return 0;let b=1,r=1;for(let i=1;i<s.length;i++){r=s[i]===s[i-1]?r+1:1;if(r>b)b=r;}return b;};
  function repUnits(s,u){let b=1;for(let i=0;i<s.length-u;i++){const m=s.substr(i,u);let c=1,j=i+u;while(s.substr(j,u)===m){c++;j+=u;}if(c>b)b=c;}return b;}
  const clamp=s=>{const t=s.slice(-5);let n=0;for(const c of t)if(c==='G'||c==='C')n++;return n;};
  function shannon(s,k){const f={};let n=0;for(let i=0;i+k<=s.length;i++){const w=s.substr(i,k);f[w]=(f[w]||0)+1;n++;}let e=0;for(const w in f){const p=f[w]/n;e-=p*Math.log2(p);}return e;}
  const g4=s=>/G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}/.test(s)||/C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}[ACGT]{1,7}C{3,}/.test(s);
  function gcWindowMax(s,w=50){if(s.length<=w)return gcf(s);let n=0;for(let i=0;i<w;i++)if('GC'.includes(s[i]))n++;let m=n;for(let i=w;i<s.length;i++){if('GC'.includes(s[i]))n++;if('GC'.includes(s[i-w]))n--;if(n>m)m=n;}return m/w;}
  // ---- design.py ----
  function enumerate(target,orient,cons,c,region,variantPenalty,onRej){
    const lo=Math.max(0,region?region[0]:0),hi=Math.min(target.length,region?region[1]:target.length);
    const acc=[];
    for(let st=lo;st<hi;st++)for(let L=cons.length_min;L<=cons.length_max;L++){
      const en=st+L;if(en>hi)break;const win=target.slice(st,en);if(/[^ACGT]/.test(win)){onRej(orient,win,st,'SEQ-04','belirsiz baz',win.match(/[^ACGT]/)[0]);continue;}
      const p=orient==='F'?win:rc(win);
      const rej=(r,why,v)=>onRej(orient,p,st,r,why,v);
      const g=gcf(p);if(!(cons.gc_min<=g&&g<=cons.gc_max)){rej('CMP-01','GC penceresi dışı',(100*g).toFixed(1)+'%');continue;}
      const cl=clamp(p);if(!(cons.gc_clamp_min<=cl&&cl<=cons.gc_clamp_max)){rej('CMP-02',"3' GC clamp",cl);continue;}
      const ho=longestHomo(p);if(ho>=cons.max_homopolymer){rej('CMP-03','homopolimer',ho);continue;}
      if(p.includes('GGGG')){rej('CMP-03','poli-G','GGGG');continue;}
      const rp=Math.max(repUnits(p,2),repUnits(p,3));if(rp>=cons.max_repeat_units){rej('CMP-04','di/tri-nükleotit tekrar',rp);continue;}
      const t=thermo(p,c,cons.salt_model);
      if(!(cons.tm_min<=t.tm&&t.tm<=cons.tm_max)){rej('TH-01','Tm penceresi dışı',t.tm.toFixed(2)+' °C');continue;}
      if(t.thetaTa<cons.theta_min){rej('TH-03','θ(Ta) eşik altında',t.thetaTa.toFixed(4));continue;}
      if(variantPenalty){const v=variantPenalty(st,en,orient);if(v){rej(v.rule,v.why,v.value);continue;}}
      acc.push({sequence:p,start:st,end:en,orientation:orient,tm_c:t.tm,gc:g,theta_ta:t.thetaTa});
    }
    acc.sort((a,b)=>Math.abs(a.tm_c-cons.tm_opt)-Math.abs(b.tm_c-cons.tm_opt));
    const out=[];const MAXD=4000;
    for(const cd of acc.slice(0,MAXD)){
      const da=bestDimer(cd.sequence,cd.sequence,'any'),d3=bestDimer(cd.sequence,cd.sequence,'3prime');
      const ga=da?da.dG:0,g3=d3?d3.dG:0;
      if(ga<cons.self_dimer_dg){onRej(orient,cd.sequence,cd.start,'DIM-01','self-dimer',ga.toFixed(2));continue;}
      if(g3<cons.self_dimer3p_dg){onRej(orient,cd.sequence,cd.start,'DIM-02',"3'-uçlu self-dimer",g3.toFixed(2));continue;}
      out.push({...cd,self_dimer_dg:ga,self_dimer3p_dg:g3});
    }
    for(const cd of acc.slice(MAXD))onRej(orient,cd.sequence,cd.start,'LIMIT','dimer taraması sınırı',MAXD);
    return out;
  }
  function buildPairs(L,Rr,cons,onRej,mustCover){
    const pairs=[];const rs=[...Rr].sort((a,b)=>a.end-b.end);
    for(const l of L){for(const r of rs){
      if(r.start<l.end)continue;const len=r.end-l.start;if(len<cons.product_min)continue;if(len>cons.product_max)break;
      if(mustCover&&!(l.start<=mustCover[0]&&r.end>=mustCover[1]))continue;
      const dtm=Math.abs(l.tm_c-r.tm_c);if(dtm>cons.max_delta_tm){onRej('P',l.sequence+'+'+r.sequence,l.start,'PAIR-01','ΔTm',dtm.toFixed(2));continue;}
      const da=bestDimer(l.sequence,r.sequence,'any'),d3=bestDimer(l.sequence,r.sequence,'3prime');const ga=da?da.dG:0,g3=d3?d3.dG:0;
      if(ga<cons.hetero_dimer_dg){onRej('P',l.sequence+'+'+r.sequence,l.start,'DIM-03','hetero-dimer',ga.toFixed(2));continue;}
      if(g3<cons.hetero_dimer3p_dg){onRej('P',l.sequence+'+'+r.sequence,l.start,'DIM-04',"3'-uçlu hetero-dimer",g3.toFixed(2));continue;}
      pairs.push({left:l,right:r,product_start:l.start,product_end:r.end,product_length:len,delta_tm:dtm,hetero_dg:ga,hetero3p_dg:g3});
      if(pairs.length>=cons.max_pairs)return pairs;
    }}
    return pairs;
  }
  function rankPairs(pairs,cons){
    for(const p of pairs){const tmd=(Math.abs(p.left.tm_c-cons.tm_opt)+Math.abs(p.right.tm_c-cons.tm_opt))/2;
      const gcd=(Math.abs(p.left.gc-.5)+Math.abs(p.right.gc-.5))/2;
      const dm=-Math.min(0,Math.min(p.hetero_dg,p.left.self_dimer_dg,p.right.self_dimer_dg));
      const d3=-Math.min(0,Math.min(p.hetero3p_dg,p.left.self_dimer3p_dg,p.right.self_dimer3p_dg));
      const th=(1-p.left.theta_ta)+(1-p.right.theta_ta);p.obj=[tmd,p.delta_tm,gcd,dm,d3,th];}
    if(pairs.length){const n=6;const lo=[],hi=[];for(let i=0;i<n;i++){lo[i]=Math.min(...pairs.map(p=>p.obj[i]));hi[i]=Math.max(...pairs.map(p=>p.obj[i]));}
      for(const p of pairs){let t=0;for(let i=0;i<n;i++){const sp=hi[i]-lo[i];t+=sp<=0?0:(p.obj[i]-lo[i])/sp;}p.score=t/n;}}
    let rem=pairs.slice();let rank=1;for(const p of pairs)p.pareto_rank=0;
    while(rem.length){const front=[];
      for(const a of rem){let dom=false;for(const b of rem){if(b===a)continue;let all=true,any=false;for(let i=0;i<6;i++){if(b.obj[i]>a.obj[i]){all=false;break;}if(b.obj[i]<a.obj[i])any=true;}if(all&&any){dom=true;break;}}if(!dom)front.push(a);}
      for(const p of front)p.pareto_rank=rank;rem=rem.filter(p=>p.pareto_rank===0);rank++;if(rank>50){for(const p of rem)p.pareto_rank=rank;break;}}
    pairs.sort((a,b)=>(a.pareto_rank-b.pareto_rank)||(a.score-b.score));return pairs;
  }
  function design(target,cons,c,opts={}){
    const rej=[];const onRej=(o,s,st,rule,why,v)=>rej.push({o,s,st,rule,why,v:String(v)});
    const L=enumerate(target,'F',cons,c,opts.leftRegion,opts.variantPenalty,onRej);
    const Rr=enumerate(target,'R',cons,c,opts.rightRegion,opts.variantPenalty,onRej);
    let pairs=buildPairs(L,Rr,cons,onRej,opts.mustCover);pairs=rankPairs(pairs,cons);
    return {target_length:target.length,left:L,right:Rr,pairs,rejections:rej};
  }
  return {setParams,rc,gcf,watsonCrick,thermo,saltSpread,dG3,halves,bestDimer,dimerList,hairpinStem,longestHomo,repUnits,clamp,shannon,g4,gcWindowMax,design,stackE,COMP,PAIR};
})();
if(typeof module!=='undefined')module.exports=OS;
