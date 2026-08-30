import h5py, numpy as np, json, base64, collections

P='/mnt/user-data/uploads/LMNA_snRNA_human_mouse_Atlas_annotated_DK_250710.h5ad'
f=h5py.File(P,'r')
obs=f['obs']

def cats(k):
    g=obs[k]
    c=[x.decode() if isinstance(x,bytes) else str(x) for x in g['categories'][:]]
    return c, g['codes'][:].astype(np.int32)

N=obs['_index'].shape[0]

# --- lineage (12, full coverage) ---
ln_c, ln_i = cats('Integrated_celltype_Human_mouse_atlas_250709')
LN_FULL={'CM':'Cardiomyocyte','FB':'Fibroblast','EC':'Endothelial','Myeloid':'Myeloid',
 'Mural':'Mural','Lymphoid':'Lymphoid','EPIC':'Epicardial','Neuro':'Neuronal',
 'LEC':'Lymphatic EC','Mast':'Mast','Proliferating':'Proliferating','Adipo':'Adipocyte'}
ln_names=np.array([LN_FULL.get(x,x) for x in ln_c])[ln_i]

# --- cell state (Celltype_integrated, full coverage) ---
ct_c, ct_i = cats('Celltype_integrated')
ct_names=np.array(ct_c)[ct_i]

# --- cohort: group (Control / LMNA DCM), full coverage ---
gp_c, gp_i = cats('group')
GP={'Control':'Control','LMNA DCM':'LMNA-DCM'}
gp_names=np.array([GP.get(x,x) for x in gp_c])[gp_i]
DS=['Control','LMNA-DCM']
ds_code=np.where(gp_names=='LMNA-DCM',1,0).astype(np.uint8)

# --- species / dataset / condition ---
sp_c, sp_i = cats('Speices'); species=np.array(sp_c)[sp_i]
dt_c, dt_i = cats('dataset'); dataset=np.array(dt_c)[dt_i]
cd_c, cd_i = cats('Condtion_integrated'); condition=np.array(cd_c)[cd_i]

# --- per-sample key (composite, best available) ---
def col(k):
    c,i=cats(k)
    a=np.array(c+[''],dtype=object)
    return a[np.where(i<0,len(c),i)]
sid=col('sample_id'); smp=col('sample'); did=col('donor_id'); oid=col('orig.ident')
sample=np.where(did!='',did,np.where(smp!='',smp,np.where(sid!='',sid,np.where(oid!='',oid,dataset))))

# --- UMAP -> uint16 ---
um=f['obsm/X_umap'][:]
x=um[:,0]; y=um[:,1]
def norm(a):
    lo,hi=np.percentile(a,0.05),np.percentile(a,99.95)
    a=np.clip(a,lo,hi)
    return ((a-a.min())/(a.max()-a.min())*65535).astype(np.uint16)
xn,yn=norm(x),norm(y)
xr=(np.clip(x,np.percentile(x,.05),np.percentile(x,99.95)))
yr=(np.clip(y,np.percentile(y,.05),np.percentile(y,99.95)))
ASPECT=float((xr.max()-xr.min())/(yr.max()-yr.min()))

xy=np.empty(N*2,dtype=np.uint16); xy[0::2]=xn; xy[1::2]=yn

# --- categorical codes as uint8 ---
ct_cats=[c for c,_ in collections.Counter(ct_names).most_common()]
ln_cats=[c for c,_ in collections.Counter(ln_names).most_common()]
ct_map={c:i for i,c in enumerate(ct_cats)}; ln_map={c:i for i,c in enumerate(ln_cats)}
assert len(ct_cats)<256 and len(ln_cats)<256
ctA=np.array([ct_map[c] for c in ct_names],dtype=np.uint8)
lnA=np.array([ln_map[c] for c in ln_names],dtype=np.uint8)

# --- gene panel ---
var=[x.decode() for x in f['var/_index'][:]]
vidx={g:i for i,g in enumerate(var)}
PANEL=['LMNA','PARP1','NAMPT','SIRT1','H2AX','TP53BP1','ATM','CDKN1A','TP53',
       'NPPA','NPPB','MYH7','TTN','RYR2','TNNT2',
       'POSTN','COL1A1','FN1','DCN','LUM','TGFB1','TIMP1','ACTA2','MYH11','PDGFRB',
       'PECAM1','CDH5','VWF','LYVE1','PTPRC','ADIPOQ','NRXN1','MKI67','MEG3','NNMT']
PANEL=[g for g in PANEL if g in vidx]
G=len(PANEL)
ip=f['X/indptr'][:]; Xd=f['X/data']; Xi=f['X/indices']
expr=np.zeros((G,N),dtype=np.uint8); gmax=[]
dense=np.zeros((G,N),dtype=np.float32)
for k,g in enumerate(PANEL):
    j=vidx[g]; a,b=ip[j],ip[j+1]
    col_v=np.zeros(N,dtype=np.float32)
    if b>a:
        col_v[Xi[a:b]]=Xd[a:b]
    dense[k]=col_v
    m=float(np.percentile(col_v[col_v>0],99.5)) if (col_v>0).any() else 1.0
    m=max(m,1e-6); gmax.append(round(m,3))
    expr[k]=np.clip(col_v/m,0,1).mul(255) if False else (np.clip(col_v/m,0,1)*255).astype(np.uint8)

# --- stats ---
def counts(arr):
    c=collections.Counter(arr)
    return [{'name':k,'n':int(v)} for k,v in c.most_common()]

comp=[]
for lname in [d['name'] for d in counts(ln_names)]:
    row={'name':lname}
    for ds in DS:
        sel=(gp_names==ds)
        row[ds]=float((ln_names[sel]==lname).mean())
    comp.append(row)

donors=[]
for s in sorted(set(sample)):
    sel=sample==s
    donors.append({'name':str(s),'n':int(sel.sum()),
        'species':str(collections.Counter(species[sel]).most_common(1)[0][0]),
        'group':str(collections.Counter(gp_names[sel]).most_common(1)[0][0]),
        'dataset':str(collections.Counter(dataset[sel]).most_common(1)[0][0])})
donors.sort(key=lambda d:-d['n'])

mean_by_type={}
for t in ct_cats:
    sel=ct_names==t
    if sel.sum()<20: continue
    mean_by_type[t]=[round(float(dense[k][sel].mean()),3) for k in range(G)]

stats={
 'n_cells':int(N),'n_genes':int(len(var)),'n_donors':len(donors),
 'n_types':len(ct_cats),'n_lineages':len(ln_cats),
 'cell_types':counts(ct_names),'lineages':counts(ln_names),
 'datasets':[{'name':ds,'n':int((gp_names==ds).sum())} for ds in DS],
 'species':counts(species),'sources':counts(dataset),'conditions':counts(condition),
 'donors':donors,'composition_lineage':comp,
 'genes':PANEL,'mean_by_type':mean_by_type,
}

blob=b''.join([xy.tobytes(),ctA.tobytes(),lnA.tobytes(),ds_code.tobytes(),expr.tobytes()])
A={'stats':stats,'n':int(N),'ct_cats':ct_cats,'ln_cats':ln_cats,'ds_cats':DS,
   'genes':PANEL,'gmax':gmax,'aspect':round(ASPECT,4),
   'blob':base64.b64encode(blob).decode()}
json.dump(A,open('atlas.json','w'),separators=(',',':'))
print('N',N,'G',G,'states',len(ct_cats),'lineages',len(ln_cats),'donors',len(donors))
print('aspect',ASPECT,'blob MB',len(blob)/1e6)
print('cohorts',stats['datasets'])
print('species',stats['species'])
print('sources',stats['sources'])
