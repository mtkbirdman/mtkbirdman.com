import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';

const css = readFileSync(new URL('../../assets/css/app.css', import.meta.url), 'utf8');
const php = readFileSync(new URL('../../glider-landing-planner.php', import.meta.url), 'utf8');

test('v0.1.8 removes transformed containing block on compact screens', () => {
  assert.match(css, /@media\s*\(max-width:760px\),\s*\(max-width:1050px\)\s*and\s*\(orientation:portrait\)/);
  assert.match(css, /\.glp-app-root\{width:calc\(100vw\s*-\s*12px\);transform:none!important;margin-left:calc\(6px\s*-\s*50vw\)!important\}/);
  assert.match(css, /\.glp-mobile-editor:not\(\[hidden\]\)\{[^}]*position:fixed;[^}]*bottom:0;/);
  assert.match(php, /Version:\s*1\.0\.0/);
});

test('v0.1.8 preserves the desktop centering and mobile controls', () => {
  assert.match(css, /\.glp-app-root\{[^}]*transform:translateX\(-50%\)/);
  assert.match(css, /\.glp-mobile-editor\[hidden\]\{display:none!important\}/);
  assert.match(css, /\.glp-mobile-editor-content\{overflow-y:auto;/);
});
