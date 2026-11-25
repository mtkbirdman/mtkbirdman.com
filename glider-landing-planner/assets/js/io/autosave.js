const DB = 'glp-autosave-v1', STORE = 'snapshots', BACKUPS = 'legacy-backups';
function openDb() {
    return new Promise((res, rej) => {
        if (!('indexedDB' in window)) {
            rej(Error('IndexedDBが利用できません。'));
            return;
        }
        const r = indexedDB.open(DB, 2);
        r.onupgradeneeded = () => { if (!r.result.objectStoreNames.contains(STORE))
            r.result.createObjectStore(STORE, { keyPath: 'savedAt' }); if (!r.result.objectStoreNames.contains(BACKUPS))
            r.result.createObjectStore(BACKUPS); };
        r.onsuccess = () => res(r.result);
        r.onerror = () => rej(r.error ?? Error('DB open failed'));
    });
}
export async function saveSnapshot(project) {
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
/** Protected copy outside the 10-snapshot rolling history, before any model-A approval. */
export async function saveLegacyBackup(project) {
    const db = await openDb();
    try {
        await new Promise((resolve, reject) => {
            const tx = db.transaction(BACKUPS, 'readwrite');
            tx.objectStore(BACKUPS).put(structuredClone(project), 'latest-model-A');
            tx.oncomplete = () => resolve();
            tx.onerror = () => reject(tx.error ?? Error('旧計画のバックアップに失敗しました。'));
        });
    }
    finally {
        db.close();
    }
}
export async function latestLegacyBackup() {
    const db = await openDb();
    try {
        return await new Promise((resolve, reject) => {
            const req = db.transaction(BACKUPS, 'readonly').objectStore(BACKUPS).get('latest-model-A');
            req.onsuccess = () => resolve(req.result ?? null);
            req.onerror = () => reject(req.error);
        });
    }
    finally {
        db.close();
    }
}
