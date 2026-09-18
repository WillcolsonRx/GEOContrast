"""Render the synthetic README animation. Requires Pillow and Node.js."""
from pathlib import Path
import subprocess,json,math
from PIL import Image,ImageDraw,ImageFont
root=Path(__file__).resolve().parents[1]
rows=json.loads(subprocess.check_output(['node','-e',"console.log(JSON.stringify(require('./docs/assets/contrast-core.js').createResults()))"],cwd=root))
a=root/'docs/assets'
font='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
def f(size):return ImageFont.truetype(font,size)
frames=[]
for frame in range(40):
 im=Image.new('RGB',(1000,500),'#101820');d=ImageDraw.Draw(im);yaw=frame/40*math.tau+.45;pitch=.28
 def proj(x,y,z):
  rx=x*math.cos(yaw)+z*math.sin(yaw);rz=-x*math.sin(yaw)+z*math.cos(yaw);ry=y*math.cos(pitch)-rz*math.sin(pitch);dep=y*math.sin(pitch)+rz*math.cos(pitch);s=7/(7-dep)
  return (610+rx*85*s,240-ry*85*s),dep,s
 d.text((30,25),'GEOContrast',font=f(23),fill='#e5e9eb');d.text((30,60),'SYNTHETIC RESULTS / 3D WEBSITE EXPLORER',font=f(11),fill='#a1b4c3')
 d.text((30,162),'A different angle.',font=f(22),fill='#c9d6ab');d.text((30,199),'The same comparison.',font=f(17),fill='#a1b4c3')
 for idx,(label,c) in enumerate([('Higher','#f0906d'),('Lower','#79a5f5'),('Not passing','#647586')]):
  y=270+idx*30;d.ellipse((31,y,39,y+8),fill=c);d.text((51,y-5),label,font=f(13),fill='#a1b4c3')
 d.text((30,407),'Adjusted P < 0.05',font=f(12),fill='#a1b4c3');d.text((30,428),'|log₂ fold change| ≥ 1',font=f(12),fill='#a1b4c3')
 for k in range(-2,3):
  d.line([proj(k,-1.5,-1.5)[0],proj(k,-1.5,1.5)[0]],fill='#2d3b48');d.line([proj(-2,-1.5,k*.75)[0],proj(2,-1.5,k*.75)[0]],fill='#2d3b48')
 points=[]
 for r in rows:
  xy,dep,s=proj(r['fc']/3,-math.log10(r['q'])/5*3-1.5,(math.log10(r['mean'])-2.925)/1.25)
  col=('#f0906d' if r['fc']>=0 else '#79a5f5') if r['q']<.05 and abs(r['fc'])>=1 else '#647586'
  points.append((dep,xy,3.3*s,col))
 for dep,(x,y),rad,col in sorted(points):d.ellipse((x-rad,y-rad,x+rad,y+rad),fill=col)
 d.text((350,467),'EFFECT × ADJUSTED SIGNIFICANCE × ABUNDANCE',font=f(11),fill='#a1b4c3')
 if frame==0:im.save(a/'landscape-still.png')
 frames.append(im.quantize(colors=96))
frames[0].save(a/'landscape.gif',save_all=True,append_images=frames[1:],duration=100,loop=0,optimize=True)
social=Image.new('RGB',(1280,640),'#101820');social.paste(Image.open(a/'landscape-still.png'),(280,135));d=ImageDraw.Draw(social);d.text((55,40),'GEOContrast',font=f(60),fill='#e5e9eb');d.text((59,117),'Make the comparison. Keep the context.',font=f(24),fill='#c9d6ab');social.save(a/'social-preview.png')
print('Rendered synthetic animation, static alternative and social card.')
