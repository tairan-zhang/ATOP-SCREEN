import {Button} from "@atop/ui";
import { useEffect, useState } from "react";
import { desktopApi } from "../services/desktop-api";
import type { Artifact } from "../types/atop";
import { Icon } from "./Icon";
import { ProcessingIndicator } from "./ProcessingIndicator";

type ArtifactPanelProps = {
  artifact: Artifact | null;
  projectName: string;
  projectPath: string;
};

function downloadTimestamp(date = new Date()) {
  const part = (value: number) => String(value).padStart(2, "0");
  return `${date.getFullYear()}_${part(date.getMonth() + 1)}_${part(date.getDate())}_${part(date.getHours())}_${part(date.getMinutes())}_${part(date.getSeconds())}`;
}

export function ArtifactPanel({ artifact, projectName, projectPath }: ArtifactPanelProps) {
  const [activeUrl, setActiveUrl] = useState("");
  const [loading, setLoading] = useState(false);
  const [downloading, setDownloading] = useState(false);
  const [error, setError] = useState("");

  useEffect(() => {
    const pngPath = artifact?.paths.png;
    if (!pngPath) {
      setActiveUrl("");
      setLoading(false);
      setError("");
      return;
    }
    let cancelled = false;
    let objectUrl = "";
    setLoading(true);
    setError("");
    setActiveUrl("");
    void desktopApi.readAssetBytes(pngPath).then(({ bytes, mimeType }) => {
      const previewBytes = new Uint8Array(bytes);
      if (cancelled) return;
      objectUrl = URL.createObjectURL(new Blob([previewBytes.buffer], { type: mimeType }));
      if (!cancelled) {
        setActiveUrl(objectUrl);
      }
    }).catch((requestError: unknown) => {
      if (!cancelled) setError(requestError instanceof Error ? requestError.message : String(requestError));
    }).finally(() => {
      if (!cancelled) setLoading(false);
    });
    return () => {
      cancelled = true;
      if (objectUrl) URL.revokeObjectURL(objectUrl);
    };
  }, [artifact]);

  async function downloadArtifactSet() {
    if (!artifact || downloading) return;
    setDownloading(true);
    setError("");
    try {
      await desktopApi.downloadArtifacts(artifact.paths, `${projectName}_${artifact.downloadLabel || artifact.title}_${downloadTimestamp()}`, projectPath);
    } catch (requestError) {
      setError(requestError instanceof Error ? requestError.message : String(requestError));
    } finally {
      setDownloading(false);
    }
  }

  return <section className="figure-output" aria-label="Figure preview">
    <header className="canvas-heading artifact-action-bar">
      <Icon name="image"/><h2>Figure preview</h2>
      {artifact && <Button variant="secondary" aria-label="Download PNG, SVG, and PDF"
        title="Download PNG, SVG and PDF" disabled={downloading} onClick={() => void downloadArtifactSet()} type="button">
        <Icon name="download" size={14}/>{downloading ? "Saving…" : "Download"}
      </Button>}
    </header>
    <div className="artifact-body">
      {!artifact ? <div className="artifact-empty"><Icon name="image" size={28}/><span>Render a figure to preview it.</span></div> :
        loading ? <ProcessingIndicator label="Loading preview…"/> :
        error ? <div role="alert" className="artifact-error">{error}</div> :
        activeUrl ? <div className="image-preview"><img alt={artifact.title} src={activeUrl}/></div> :
        <div className="artifact-empty">PNG preview unavailable</div>}
    </div>
  </section>;
}
