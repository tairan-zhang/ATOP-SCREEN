const artwork: Record<string, string> = {
  "activity": "ChartNoAxesCombined",
  "arrowLeft": "ArrowLeft",
  "arrowRight": "ArrowRight",
  "arrowUp": "ArrowUp",
  "bar": "ChartNoAxesCombined",
  "check": "Check",
  "chevronDown": "ChevronDown",
  "chevronRight": "ChevronRight",
  "close": "X",
  "database": "Sheet",
  "download": "Download",
  "file": "Files",
  "flask": "Dna",
  "folder": "Folder",
  "image": "Image",
  "info": "Info",
  "moon": "PanelsTopLeft",
  "panel": "PanelLeft",
  "panelRight": "PanelRight",
  "pdf": "FileText",
  "play": "ArrowRight",
  "plus": "Plus",
  "refresh": "RefreshCw",
  "settings": "Settings2",
  "sun": "PanelsTopLeft",
  "upload": "Upload",
  "trash": "Archive",
  "home": "House",
  "newAnalysis": "SquarePen"
};

export function Icon({ name, size = 16 }: { name: string; size?: number }) {
  return <img src={import.meta.env.BASE_URL + "ui-icons/" + (artwork[name] ?? "Info") + ".svg"}
    alt="" aria-hidden="true" draggable={false} className="library-icon" style={{width:size,height:size}} />;
}
