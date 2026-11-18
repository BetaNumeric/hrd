const appRoot = new URL(self.registration.scope);
const cachePrefix = `hr-diagram:${appRoot.href}:`;
const cacheName = `${cachePrefix}v1`;
const appShell = new URL("index.html", appRoot);
const assets = ["manifest.json", "icons/icon_32.png", "icons/icon_192.png", "icons/icon_512.png"]
  .map((path) => new URL(path, appRoot));
const cacheKeys = new Map(assets.map((url) => [url.pathname, url.href]));
cacheKeys.set(appRoot.pathname, appShell.href);
cacheKeys.set(appShell.pathname, appShell.href);

self.addEventListener("install", (event) => {
  // The single HTML file already contains the catalogs, fonts, and renderer.
  // Store one copy for both the launch URL and direct index.html visits.
  const requests = [appShell, ...assets].map((url) => new Request(url, { cache: "no-cache" }));
  event.waitUntil((async () => {
    const cache = await caches.open(cacheName);
    await cache.addAll(requests);
    await self.skipWaiting();
  })());
});

self.addEventListener("activate", (event) => {
  event.waitUntil((async () => {
    const names = await caches.keys();
    await Promise.all(names.filter((name) => name.startsWith(cachePrefix) && name !== cacheName)
      .map((name) => caches.delete(name)));
    await self.clients.claim();
  })());
});

self.addEventListener("fetch", (event) => {
  const { request } = event;
  if (request.method !== "GET") return;
  const url = new URL(request.url);
  if (url.origin !== appRoot.origin) return;
  const key = cacheKeys.get(url.pathname);
  if (!key) return;

  // Revalidate online so ordinary HTML and catalog updates arrive on reload.
  // Deliver the response while the cache write finishes in the background.
  const network = fetch(request, { cache: "no-cache" }).catch(() => null);
  event.waitUntil(network.then(async (response) => {
    if (!response?.ok) return;
    const copy = response.clone();
    const cache = await caches.open(cacheName);
    await cache.put(key, copy);
  }).catch(() => {}));
  event.respondWith(network.then(async (response) => {
    if (response?.ok) return response;
    const cached = await caches.match(key, { cacheName }).catch(() => null);
    return cached ?? response ?? Response.error();
  }));
});
