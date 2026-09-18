/** Dependency-free local preview. R Shiny runs separately. */
import http from 'node:http';
import path from 'node:path';
import fs from 'node:fs/promises';
import {fileURLToPath} from 'node:url';
const root=fileURLToPath(new URL('../docs/',import.meta.url));
const port=Number(process.env.PORT||4173);
const types={'.html':'text/html; charset=utf-8','.js':'text/javascript; charset=utf-8','.css':'text/css; charset=utf-8','.svg':'image/svg+xml','.png':'image/png','.gif':'image/gif','.txt':'text/plain; charset=utf-8','.zip':'application/zip'};
const server=http.createServer(async(req,res)=>{
 if(!['GET','HEAD'].includes(req.method)){res.writeHead(405,{'Allow':'GET, HEAD'});res.end();return;}
 try{
  const url=new URL(req.url,'http://localhost');
  const requested=decodeURIComponent(url.pathname);
  const file=path.resolve(root,'.'+(requested.endsWith('/')?requested+'index.html':requested));
  if(!file.startsWith(root)){res.writeHead(403);res.end('Forbidden');return;}
  const stat=await fs.stat(file);if(!stat.isFile())throw Error('Not a file');
  const body=req.method==='HEAD'?null:await fs.readFile(file);
  res.writeHead(200,{'Content-Type':types[path.extname(file)]||'application/octet-stream','Content-Length':stat.size,'X-Content-Type-Options':'nosniff','Cache-Control':'no-cache'});res.end(body);
 }catch{res.writeHead(404,{'Content-Type':'text/plain; charset=utf-8'});res.end('Not found');}
});
server.on('error',error=>{console.error(error.message);process.exitCode=1;});
server.listen(port,'127.0.0.1',()=>console.log(`GEOContrast preview: http://127.0.0.1:${server.address().port}`));
