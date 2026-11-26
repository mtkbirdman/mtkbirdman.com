import { validateProject } from './project.js?v=0.1.20';
const DB = 'glp-autosave-v1', STORE = 'snapshots';
function openDb() {
    return new Promise((res, rej) => {
        if (!('indexedDB' in window)) {
            rej(Error('IndexedDBが利用できません。'));
            return;
        }
        const r = indexedDB.open(DB, 2);
        r.onupgradeneeded = () => { if (!r.result.objectStoreNames.contains(STORE))
            r.result.createObjectStore(STORE, { keyPath: 'savedAt' }); };
        r.onsuccess = () => res(r.result);
        r.onerror = () => rej(r.error ?? Error('DB open failed'));
    });
}
export async function saveSnapshot(project) {
    validateProject(project);
    const db = await openDb();
    try {
        await new Promise((res, rej) => { const tx = db.transaction(STORE, 'readwrite'), store = tx.objectStore(STORE); store.put({ savedAt: Date.now(), project }); tx.oncomplete = () => res(); tx.onerror = () => rej(tx.error ?? Error('DB save failed')); });
        await new Promise((res, rej) => { const tx = db.transaction(STORE, 'readwrite'), store = tx.objectStore(STORE), req = store.getAllKeys(); req.onsuccess = () => { const keys = req.result; if (keys.length > 10)
            keys.slice(0, keys.length - 10).forEach(k => store.delete(k)); }; tx.oncomplete = () => res(); tx.onerror = () => rej(tx.error ?? Error('DB purge failed')); });
    }
    finally {
        db.close();
    }
}
export async function latestSnapshot() { const db = await openDb(); try {
    return await new Promise((res, rej) => { const tx = db.transaction(STORE, 'readonly'), req = tx.objectStore(STORE).openCursor(null, 'prev'); req.onsuccess = () => res(req.result ? req.result.value.project : null); req.onerror = () => rej(req.error); });
}
finally {
    db.close();
} }
