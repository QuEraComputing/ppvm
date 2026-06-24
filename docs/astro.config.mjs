import { defineConfig } from "astro/config";

// `site` + `base` are read from env in CI so the same astro config powers
// three deploy targets. The repo is public, so Pages serves from the
// `queracomputing.github.io/ppvm/` project page — every deploy except
// local dev is therefore based under `/ppvm/`:
//   · local dev:           PPVM_SITE=http://localhost:4321              PPVM_BASE=/
//   · main → gh-pages:     PPVM_SITE=https://queracomputing.github.io   PPVM_BASE=/ppvm/
//   · PR preview deploy:   PPVM_SITE=https://queracomputing.github.io   PPVM_BASE=/ppvm/pr-preview/pr-<N>
const site = process.env.PPVM_SITE ?? "https://queracomputing.github.io";
const base = process.env.PPVM_BASE ?? "/";

export default defineConfig({
  site,
  base,
  trailingSlash: "ignore",
  devToolbar: { enabled: false },
});
