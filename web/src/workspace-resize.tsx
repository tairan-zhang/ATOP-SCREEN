import React,{useEffect,useRef,type PointerEvent as ReactPointerEvent} from 'react';
import {flushSync} from 'react-dom';
import {Button} from '@atop/ui';

/** Resize the actual panel each frame; persist preferences only at gesture end. */
export function WorkspaceResizeHandle({label,controls,className='',target,variable,edge,limits,onCommit,onReset}:{
 label:string;controls:string;className?:string;target:()=>HTMLElement|null;variable:string;edge:'left'|'right';
 limits:()=>{min:number;max:number};onCommit:(width:number)=>void;onReset:()=>void;
}){
 const callbacks=useRef({target,limits,onCommit,onReset});callbacks.current={target,limits,onCommit,onReset};
 const stop=useRef<((commit:boolean)=>void)|undefined>(undefined);const handle=useRef<HTMLButtonElement>(null);
 const sign=edge==='right'?1:-1;
 function clamp(value:number,bounds=callbacks.current.limits()){return Math.round(Math.min(bounds.max,Math.max(bounds.min,value)));}
 function announce(value:number){handle.current?.setAttribute('aria-valuenow',String(Math.round(value)));handle.current?.setAttribute('aria-valuetext',Math.round(value)+' 像素');}
 useEffect(()=>{
  const panel=callbacks.current.target();if(!panel)return;
  // Observer dimensions are already computed by the browser: no layout read/state update.
  const observer=new ResizeObserver(entries=>{if(stop.current)return;const box=entries[0]?.borderBoxSize?.[0];if(box)announce(box.inlineSize);});
  observer.observe(panel);const bounds=callbacks.current.limits();handle.current?.setAttribute('aria-valuemin',String(bounds.min));handle.current?.setAttribute('aria-valuemax',String(bounds.max));announce(panel.getBoundingClientRect().width);
  return()=>{observer.disconnect();stop.current?.(false);};
 },[]);
 function start(event:ReactPointerEvent<HTMLButtonElement>){
  if(!event.pointerType||event.button!==0||stop.current||document.documentElement.dataset.atopResizing)return;
  const panel=callbacks.current.target();if(!panel)return;
  const button=event.currentTarget,pointer=event.pointerId,startX=event.clientX,initial=panel.getBoundingClientRect().width;
  const bounds=callbacks.current.limits();let next=initial,frame=0,ended=false,moved=false;
  const oldValue=panel.style.getPropertyValue(variable),oldPriority=panel.style.getPropertyPriority(variable);
  const oldCursor=document.documentElement.style.cursor;
  event.preventDefault();button.focus({preventScroll:true});
  try{button.setPointerCapture(pointer);}catch{/* Window listeners cover capture failures. */}
  document.documentElement.dataset.atopResizing='true';document.documentElement.style.cursor='col-resize';
  panel.dataset.liveResize='true';button.classList.add('is-resizing');
  button.setAttribute('aria-valuemin',String(bounds.min));button.setAttribute('aria-valuemax',String(bounds.max));
  window.dispatchEvent(new Event('atop-resize-state'));
  const apply=()=>{frame=0;panel.style.setProperty(variable,next+'px');announce(next);};
  const finish=(commit:boolean)=>{
   if(ended)return;ended=true;cancelAnimationFrame(frame);
   window.removeEventListener('pointermove',move);window.removeEventListener('pointerup',up);window.removeEventListener('pointercancel',cancel);window.removeEventListener('blur',blur);window.removeEventListener('keydown',key,true);window.removeEventListener('resize',cancelResize);button.removeEventListener('lostpointercapture',lost);
   try{button.releasePointerCapture(pointer);}catch{/* Capture may already be released. */}
   const value=commit&&moved?clamp(next):Math.round(initial);
   if(commit&&moved){apply();flushSync(()=>callbacks.current.onCommit(value));}
   // The saved React/CSS value now owns width. Restore only our temporary override.
   if(oldValue)panel.style.setProperty(variable,oldValue,oldPriority);else panel.style.removeProperty(variable);
   stop.current=undefined;announce(value);button.classList.remove('is-resizing');delete panel.dataset.liveResize;
   delete document.documentElement.dataset.atopResizing;document.documentElement.style.cursor=oldCursor;
   window.dispatchEvent(new Event('atop-resize-state'));
  };
  const position=(x:number)=>{const value=clamp(initial+sign*(x-startX),bounds);if(x!==startX)moved=true;if(value===next)return;next=value;if(!frame)frame=requestAnimationFrame(apply);};
  const move=(e:PointerEvent)=>{if(e.pointerId===pointer)position(e.clientX);};
  const up=(e:PointerEvent)=>{if(e.pointerId!==pointer)return;position(e.clientX);finish(true);};
  const cancel=(e:PointerEvent)=>{if(e.pointerId===pointer)finish(false);};
  const lost=()=>finish(false);const blur=()=>finish(true);const cancelResize=()=>finish(false);
  const key=(e:KeyboardEvent)=>{if(e.key==='Escape'){e.preventDefault();e.stopPropagation();finish(false);}};
  stop.current=finish;
  window.addEventListener('pointermove',move);window.addEventListener('pointerup',up);window.addEventListener('pointercancel',cancel);window.addEventListener('blur',blur);window.addEventListener('resize',cancelResize);window.addEventListener('keydown',key,true);button.addEventListener('lostpointercapture',lost);
 }
 function key(event:React.KeyboardEvent<HTMLButtonElement>){
  const width=callbacks.current.target()?.getBoundingClientRect().width;if(width===undefined||stop.current)return;
  const bounds=callbacks.current.limits();const step=event.shiftKey?32:8;
  const next=event.key==='ArrowLeft'?width-sign*step:event.key==='ArrowRight'?width+sign*step:event.key==='Home'?bounds.min:event.key==='End'?bounds.max:undefined;
  if(next!==undefined){event.preventDefault();event.stopPropagation();callbacks.current.onCommit(clamp(next,bounds));}
  else if(event.key==='Enter'){event.preventDefault();event.stopPropagation();callbacks.current.onReset();}
 }
 return <Button ref={handle} variant="ghost" pressScale={1} role="separator" aria-label={label} aria-controls={controls} aria-orientation="vertical" title={label+'；方向键微调，双击或 Enter 恢复默认'} className={'workspace-resize-handle '+className} onPointerDown={start} onKeyDownCapture={key} onDoubleClick={()=>callbacks.current.onReset()}/>;
}
