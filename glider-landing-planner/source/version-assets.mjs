// Cache-bust every relative ES module, not just the WordPress entry script.
// Keep the TypeScript sources unversioned; this is the final step after tsc.
import {readFileSync, writeFileSync, readdirSync} from 'node:fs';
import {resolve, dirname, join} from 'node:path';
import {fileURLToPath} from 'node:url';
const sourceDir = dirname(fileURLToPath(import.meta.url));
const assets = resolve(sourceDir, '../assets/js');
const version = JSON.parse(readFileSync(join(sourceDir, 'package.json'), 'utf8')).version;
let count = 0;
function visit(dir) {
    for (const entry of readdirSync(dir, {withFileTypes: true})) {
        const path = join(dir, entry.name);
        if (entry.isDirectory()) { visit(path); continue; }
        if (!entry.name.endsWith('.js')) continue;
        const original = readFileSync(path, 'utf8');
        const updated = original.replace(/(\bfrom\s*['"])(\.{1,2}\/[^'"]+\.js)(\?v=[^'"]+)?(['"])/g,
            (_, before, module, _oldVersion, after) => `${before}${module}?v=${version}${after}`);
        writeFileSync(path, updated, 'utf8');
        count++;
    }
}
visit(assets);
console.log(`Versioned relative module imports in ${count} compiled JS files (v${version}).`);
