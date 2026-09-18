const {test}=require('node:test');
const assert=require('node:assert/strict');
const {spawn}=require('node:child_process');
const path=require('node:path');
test('Node server delivers the site, assets and downloads with correct types',{timeout:10000},async()=>{
 const child=spawn(process.execPath,[path.join(__dirname,'../tools/serve.mjs')],{env:{...process.env,PORT:'0'},stdio:['ignore','pipe','pipe']});
 try{
  const base=await new Promise((resolve,reject)=>{const timer=setTimeout(()=>reject(Error('Server startup timed out')),4000);child.once('error',reject);child.stderr.on('data',data=>{clearTimeout(timer);reject(Error(String(data)));});child.stdout.on('data',data=>{const match=String(data).match(/http:\/\/127\.0\.0\.1:\d+/);if(match){clearTimeout(timer);resolve(match[0]);}});});
  for(const [url,status,mime]of[['/',200,'text/html'],['/guide.html',200,'text/html'],['/assets/contrast.js',200,'text/javascript'],['/assets/contrast-core.js',200,'text/javascript'],['/downloads/GEOContrast-app.zip',200,'application/zip'],['/missing-file',404,'text/plain']]){
   const response=await fetch(base+url);assert.equal(response.status,status,url);assert.ok(response.headers.get('content-type').includes(mime));await response.arrayBuffer();
  }
  const blocked=await fetch(base+'/%2f..%2f..%2fapp.R');assert.equal(blocked.status,403);await blocked.text();
 }finally{child.kill();}
});
