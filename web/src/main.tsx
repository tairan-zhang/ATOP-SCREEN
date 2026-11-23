import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import App from "./app/App";
import "@atop/ui/styles.css";
import "./styles/layout.css";
import "./styles/product-shell.css";
import "./styles/interaction-shapes.css";
import "./styles/screen.css";

const platform = navigator.platform.toLowerCase();
document.documentElement.classList.add(
  platform.includes("mac")
    ? "desktop-macos"
    : platform.includes("win")
      ? "platform-win"
      : platform.includes("linux")
        ? "platform-linux"
        : "platform-other",
);

document.addEventListener("contextmenu", (event) => event.preventDefault());

const root = document.getElementById("root");
if (!root) throw new Error("ATOP SCREEN root element is missing");

createRoot(root).render(
  <StrictMode>
    <App />
  </StrictMode>,
);
