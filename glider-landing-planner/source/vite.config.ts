// Optional future development build. The distributed test ZIP is compiled by `tsc`.
// `npm run build:vite` requires obtaining Vite from npm first.
import {defineConfig} from 'vite';
export default defineConfig({
  base:'./',
  build:{
    target:'es2022',
    outDir:'../assets',
    emptyOutDir:false,
    rollupOptions:{
      input:'src/main.ts',
      output:{format:'es',entryFileNames:'js/main.js',chunkFileNames:'js/[name]-[hash].js'}
    }
  }
});
