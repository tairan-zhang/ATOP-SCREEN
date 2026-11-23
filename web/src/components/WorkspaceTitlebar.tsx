import {Button} from "@atop/ui";
import {Icon} from "./Icon";

type Props = {
  projectName: string; pageTitle: string; sidebarOpen: boolean; previewOpen: boolean;
  onToggleSidebar: () => void; onTogglePreview: () => void;
};
export function WorkspaceTitlebar(props: Props) {
  return <header className="product-titlebar">
    <div className="product-window-actions">
      <Button variant="ghost" size="icon" aria-label={props.sidebarOpen ? "Hide navigation" : "Show navigation"} aria-expanded={props.sidebarOpen} onClick={props.onToggleSidebar}><Icon name="panel"/></Button>
    </div>
    <div className="screen-title"><span>{props.projectName}</span><Icon name="chevronRight" size={13}/><h1>{props.pageTitle}</h1></div>
    <div className="product-window-actions">
      <Button variant="ghost" size="icon" aria-label={props.previewOpen ? "Hide preview" : "Show preview"} aria-expanded={props.previewOpen} onClick={props.onTogglePreview}><Icon name="panelRight"/></Button>
    </div>
  </header>;
}
