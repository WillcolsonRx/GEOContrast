from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
root=Path(__file__).resolve().parents[1]
files=[root/n for n in ('app.R','install_packages.R','install_annotation_dbs.R','GEOContrast.Rproj','LICENSE')]
for folder in ('R','www','tests'):files+=list((root/folder).rglob('*'))
with ZipFile(root/'docs/downloads/GEOContrast-app.zip','w',ZIP_DEFLATED) as z:
 for p in files:
  if p.is_file():z.write(p,'GEOContrast/'+str(p.relative_to(root)))
 z.writestr('GEOContrast/README.md','# GEOContrast\n\nOpen GEOContrast.Rproj in RStudio. From this folder run:\n\n```r\nsource("install_packages.R")\nshiny::runApp()\n```\n\nOptional: source("install_annotation_dbs.R"). Internet is needed for installation and GEO downloads. One organism per analysis. RNA-seq requires integer raw counts, not TPM/FPKM/CPM. This download preserves the uploaded v0.4.8 analysis source; R runtime validation is still required. See the included field guide.\n')
 for name in ('guide.html','license.txt','assets/style.css','assets/mark.svg'):
  p=root/'docs'/name;data=p.read_bytes()
  if name=='guide.html':data=data.replace(b'index.html#explorer',b'#demo').replace(b'index.html',b'#main').replace(b'href="downloads/GEOContrast-app.zip" download',b'href="../README.md"')
  z.writestr('GEOContrast/guide/'+name,data)
print('App archive created.')
