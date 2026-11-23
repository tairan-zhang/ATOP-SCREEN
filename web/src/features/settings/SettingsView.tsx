import {Button} from "@atop/ui";
import type { HealthResult } from "@atop-screen/ipc-contracts";
import type { Theme } from "../../types/atop";
import { Icon } from "../../components/Icon";

export function SettingsView({ health, onRefresh, setTheme, theme }: { health: HealthResult | null; onRefresh: () => void; setTheme: (theme: Theme) => void; theme: Theme }) {
  return <div className="workspace-view"><section className="panel-card settings-card"><div className="result-summary-top"><div><h2>{!health ? "R runtime" : health.status === "ready" ? "R connected" : "R disconnected"}</h2></div><Button variant="secondary" className="secondary-button" onClick={onRefresh} type="button"><Icon name="refresh" />Refresh</Button></div><div className="settings-grid"><SettingRow label="Desktop shell" value={window.atop?.isDesktop ? `Electron · ${window.atop.platform}` : "Browser preview"} /><SettingRow label="R version" value={health?.rVersion ?? "Not checked"} /><SettingRow label="Loaded modules" value={health?.loadedModules?.join(", ") || "Not loaded"} /><SettingRow label="Analysis functions" value={health ? Object.entries(health.functions).filter(([, available]) => available).map(([name]) => name).join(", ") || "Unavailable" : "Not checked"} /></div></section><section className="panel-card appearance-card"><div><h2>Theme</h2></div><div className="theme-switcher"><Button variant="secondary" aria-pressed={theme === "dark"} className={theme === "dark" ? "theme-choice theme-choice-active" : "theme-choice"} onClick={() => setTheme("dark")} type="button"><Icon name="moon" />Dark</Button><Button variant="secondary" aria-pressed={theme === "light"} className={theme === "light" ? "theme-choice theme-choice-active" : "theme-choice"} onClick={() => setTheme("light")} type="button"><Icon name="sun" />Light</Button></div></section></div>;
}

function SettingRow({ label, value }: { label: string; value: string }) {
  return <div className="setting-row"><span>{label}</span><strong title={value}>{value}</strong></div>;
}
