import {useRef, useState, type CSSProperties, type ReactNode} from "react";
import {WorkspaceResizeHandle} from "../workspace-resize";

export function ParameterPane({children}: {children: ReactNode}) {
  const pane = useRef<HTMLElement>(null);
  const [width, setWidth] = useState(370);
  return <aside id="parameter-pane" ref={pane} className="parameter-pane" aria-label="Parameters"
    style={{"--parameter-width": width + "px"} as CSSProperties}>
    <div className="parameter-scroll screen-page">{children}</div>
    <WorkspaceResizeHandle label="Resize parameters" controls="parameter-pane" target={() => pane.current}
      edge="right" variable="--parameter-width" className="workspace-resize-sidebar"
      limits={() => ({min:300, max:Math.min(540, window.innerWidth * 0.45)})}
      onCommit={setWidth} onReset={() => setWidth(370)}/>
  </aside>;
}
