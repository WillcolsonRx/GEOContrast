/** Assemble an offline, one-file preview without dependencies. */
import fs from 'node:fs/promises';
import {fileURLToPath} from 'node:url';
const docs=new URL('../docs/',import.meta.url);
const read=name=>fs.readFile(new URL(name,docs),'utf8');
let home=await read('index.html');
const guide=await read('guide.html');
let guideBody=guide.match(/<main id="main" class="wrap">([\s\S]*?)<\/main>/)[1];
for(const[,id]of guideBody.matchAll(/id="([^"]+)"/g)){
 guideBody=guideBody.replaceAll(`id="${id}"`,`id="guide-${id}"`).replaceAll(`href="#${id}"`,`href="#guide-${id}"`).replaceAll(`data-copy="${id}"`,`data-copy="guide-${id}"`);
}
home=home.replace('</main>',`<section id="documentation" class="wrap preview-guide">${guideBody}</section></main>`);
home=home.replaceAll(/href="guide\.html#([^"]+)"/g,'href="#guide-$1"').replaceAll('href="guide.html"','href="#documentation"').replaceAll('href="index.html"','href="#explorer"').replaceAll('href="index.html#explorer"','href="#explorer"').replaceAll('href="index.html#workflow"','href="#workflow"');
const css=await read('assets/style.css');
home=home.replace('<link rel="stylesheet" href="assets/style.css">',`<style>${css}</style>`);
let scripts='';
for(const name of ['config.js','contrast-core.js','contrast.js']){
 home=home.replace(`<script defer src="assets/${name}"></script>`,'');scripts+='\n'+await read(`assets/${name}`);
}
home=home.replace('</body>',`<script>${scripts}</script></body>`);
for(const[name,mime]of[['mark.svg','image/svg+xml'],['social-preview.png','image/png']]){
 const data=await fs.readFile(new URL(`assets/${name}`,docs));home=home.replaceAll(`assets/${name}`,`data:${mime};base64,${data.toString('base64')}`);
}
const license=await read('license.txt');
home=home.replaceAll('href="license.txt"',`href="data:text/plain;charset=utf-8,${encodeURIComponent(license)}" download="LICENSE.txt"`);
const source=await fs.readFile(new URL('downloads/GEOContrast-app.zip',docs));
home=home.replaceAll('href="downloads/GEOContrast-app.zip" download',`href="data:application/zip;base64,${source.toString('base64')}" download="GEOContrast-app.zip"`);
const output=new URL('../GEOContrast-preview.html',import.meta.url);await fs.writeFile(output,home);
console.log(`Created ${fileURLToPath(output)} (${Buffer.byteLength(home).toLocaleString()} bytes)`);
