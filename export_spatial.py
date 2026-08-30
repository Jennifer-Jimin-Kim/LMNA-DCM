import h5py, numpy as np, json, glob, os, csv, base64, shutil

SAMPLES=[('Control-0','Control','V12U21-403-A1'),('Control-1','Control','V12U21-403-B1'),
         ('LMNA-1','LMNA-DCM','V12U21-403-C1'),('LMNA-2','LMNA-DCM','V12U21-403-D1')]
PANEL=['LMNA','POSTN','COL1A1','NPPA','NPPB','MYH7','TTN','ACTA2','PECAM1','DCN',
       'PARP1','NAMPT','CDKN1A','TIMP1','FN1','TNNT2','PTPRC','MYH11','VWF','NNMT']
OUT='/home/claude/site'
os.makedirs(OUT+'/img',exist_ok=True); os.makedirs(OUT+'/data',exist_ok=True)

def outs(s):
    g=glob.glob(f'/home/claude/vis/{s}/**/outs',recursive=True)
    return g[0]

res={'panel':PANEL,'samples':[]}
metrics=[]
for name,cohort,slide in SAMPLES:
    o=outs(name)
    sf=json.load(open(o+'/spatial/scalefactors_json.json'))
    pos={}
    with open(o+'/spatial/tissue_positions.csv') as fh:
        for r in csv.DictReader(fh):
            if int(r['in_tissue'])==1:
                pos[r['barcode']]=(float(r['pxl_col_in_fullres']),float(r['pxl_row_in_fullres']))
    h=h5py.File(o+'/filtered_feature_bc_matrix.h5','r')
    grp=h['matrix']
    bcs=[b.decode() for b in grp['barcodes'][:]]
    genes=[b.decode() for b in grp['features/name'][:]]
    gidx={}
    for i,g in enumerate(genes):
        gidx.setdefault(g,i)
    data=grp['data'][:]; indices=grp['indices'][:]; indptr=grp['indptr'][:]
    nspot=len(bcs)
    tot=np.zeros(nspot)
    vals={g:np.zeros(nspot) for g in PANEL}
    want={gidx[g]:g for g in PANEL if g in gidx}
    for j in range(nspot):
        a,b=indptr[j],indptr[j+1]
        idx=indices[a:b]; d=data[a:b]
        tot[j]=d.sum()
        hit=np.isin(idx,list(want.keys()))
        for ii,dd in zip(idx[hit],d[hit]):
            vals[want[ii]][j]=dd
    cf=1e4/np.maximum(tot,1)
    keep=[j for j,b in enumerate(bcs) if b in pos]
    xs=np.array([pos[bcs[j]][0] for j in keep])*sf['tissue_lowres_scalef']
    ys=np.array([pos[bcs[j]][1] for j in keep])*sf['tissue_lowres_scalef']
    ent={'name':name,'cohort':cohort,'slide':slide,
         'img':f'img/{name}.png','spot_r':sf['spot_diameter_fullres']*sf['tissue_lowres_scalef']/2,
         'x':[round(float(v),1) for v in xs],'y':[round(float(v),1) for v in ys],
         'umi':[int(tot[j]) for j in keep],'genes':{}}
    for g in PANEL:
        v=np.log1p(vals[g][keep]*cf[keep])
        m=float(np.percentile(v,99.5)) if v.max()>0 else 1.0
        m=max(m,1e-6)
        ent['genes'][g]={'max':round(m,3),
            'v':[int(x) for x in np.clip(v/m,0,1)*255]}
    res['samples'].append(ent)
    shutil.copy(o+'/spatial/tissue_lowres_image.png',f'{OUT}/img/{name}.png')
    mrow=list(csv.DictReader(open(o+'/metrics_summary.csv')))[0]
    metrics.append({'name':name,'cohort':cohort,'slide':slide,
        'spots':int(mrow['Number of Spots Under Tissue']),
        'median_genes':float(mrow['Median Genes per Spot']),
        'median_umi':float(mrow['Median UMI Counts per Spot']),
        'genes_detected':int(mrow['Total Genes Detected']),
        'mean_reads':float(mrow['Mean Reads Under Tissue per Spot']),
        'saturation':float(mrow['Sequencing Saturation'])})
    print(name,len(keep),'spots')

res['metrics']=metrics
json.dump(res,open(OUT+'/data/spatial.json','w'),separators=(',',':'))
print('spatial.json',os.path.getsize(OUT+'/data/spatial.json')/1e6,'MB')
