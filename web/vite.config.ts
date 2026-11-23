import { defineConfig } from "vite";
import { fileURLToPath } from "node:url";

export default defineConfig({
  base: "./",
  resolve: { alias: {
    "react/jsx-dev-runtime": fileURLToPath(new URL("./src/jsx-dev-compat.ts", import.meta.url)),
    "@": fileURLToPath(new URL("./src", import.meta.url)),
  } },
  server: { host: "127.0.0.1", port: 5173, strictPort: true },
  build: { outDir: "../dist", emptyOutDir: true, target: "chrome120" },
});
