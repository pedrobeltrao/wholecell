# origin: yeast yeast_cell_dev_bundle/web/test_viewer.py @ e3b7448 (no licence yet; local only, see wholecell/NOTICE.md)
"""Headless Chromium (SwiftShader WebGPU) test harness for the yeast cell viewer.
usage: python test_viewer.py <mode> [args]
modes: boot | bookmarks [i,j,...] | eval "<js>"
"""
import asyncio, sys, subprocess, time, json, os, socket
from playwright.async_api import async_playwright

ROOT = os.environ.get('YC_ROOT', '/workspace/yc3d/web/app/dist')
PORT = int(os.environ.get('YC_PORT', '8811'))
OUT = os.environ.get('YC_OUT', '/workspace/yc3d/web/shots')
W, H = int(os.environ.get('YC_W', '1280')), int(os.environ.get('YC_H', '800'))
ARGS = ["--enable-unsafe-webgpu", "--enable-unsafe-swiftshader", "--enable-gpu", "--use-gl=angle",
        "--use-angle=swiftshader", "--enable-features=Vulkan,WebGPUService", "--use-vulkan=swiftshader"]


def port_open(p):
    s = socket.socket(); r = s.connect_ex(('127.0.0.1', p)); s.close(); return r == 0


async def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'boot'
    os.makedirs(OUT, exist_ok=True)
    srv = None
    if not port_open(PORT):
        srv = subprocess.Popen([sys.executable, '-m', 'http.server', str(PORT), '-d', ROOT],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for _ in range(50):
            if port_open(PORT): break
            time.sleep(0.1)
    logs = []
    try:
        async with async_playwright() as p:
            b = await p.chromium.launch(headless=True, args=ARGS)
            pg = await b.new_page(viewport={'width': W, 'height': H})
            pg.on('console', lambda m: logs.append(f'[{m.type}] {m.text}'))
            pg.on('pageerror', lambda e: logs.append(f'[pageerror] {e}'))
            t0 = time.time()
            await pg.goto(f"http://127.0.0.1:{PORT}/{os.environ.get('YC_PAGE', 'index.html')}{os.environ.get('YC_QS', '')}")
            try:
                await pg.wait_for_function('window.yc && window.yc.ready', timeout=120000)
            except Exception as e:
                logs.append(f'[harness] not ready: {e}')
            print(f'ready after {time.time()-t0:.1f}s')
            if os.environ.get('YC_Q'):
                await pg.evaluate(f"window.yc.setQuality('{os.environ['YC_Q']}')")

            async def settle(maxs=90):
                t = time.time()
                while time.time() - t < maxs:
                    st = await pg.evaluate('window.yc && window.yc.state()')
                    if st and st.get('pending', 1) == 0 and st.get('detailPending', 0) == 0:
                        break
                    await pg.wait_for_timeout(1000)
                for _ in range(2):   # two full frames after streaming finished so stats + image are current
                    await pg.evaluate('window.yc && window.yc.render()')
                    await pg.wait_for_timeout(300)
                    t2 = time.time()
                    while time.time() - t2 < 120:
                        st2 = await pg.evaluate('window.yc.state()')
                        if not st2.get('gpuBusy') and not st2.get('dirty'): break
                        await pg.wait_for_timeout(300)
                return await pg.evaluate('window.yc && window.yc.state()')

            if mode == 'boot':
                st = await settle(int(os.environ.get('YC_SETTLE', '90')))
                print('state', json.dumps(st)[:1500])
                await pg.screenshot(path=f'{OUT}/boot.png', timeout=int(os.environ.get('YC_SHOT_TO','240000')))
            elif mode == 'bookmarks':
                sel = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else None
                n = await pg.evaluate('window.yc.bookmarks.length')
                for i in (sel if sel is not None else range(n)):
                    name = await pg.evaluate(f'window.yc.bookmarks[{i}].name')
                    await pg.evaluate(f'window.yc.applyBookmark({i}, false)')
                    t = time.time()
                    st = await settle(int(os.environ.get('YC_SETTLE', '120')))
                    fn = f"{OUT}/{i:02d}_{name.replace(' ', '_').replace('/', '-')}.png"
                    await pg.screenshot(path=fn, timeout=int(os.environ.get('YC_SHOT_TO','240000')))
                    print(i, name, f'{time.time()-t:.1f}s', json.dumps(st)[:600])
            elif mode == 'eval':
                await settle(int(os.environ.get('YC_SETTLE', '60')))
                for js in sys.argv[2:]:
                    r = await pg.evaluate(js)
                    print('EVAL', json.dumps(r)[:3000] if not isinstance(r, str) else r[:3000])
                    await pg.wait_for_timeout(500)
                await settle(int(os.environ.get('YC_SETTLE', '60')))
                await pg.screenshot(path=f'{OUT}/eval.png', timeout=int(os.environ.get('YC_SHOT_TO','240000')))
            await b.close()
    finally:
        if srv: srv.terminate()
        print('---- console (%d) ----' % len(logs))
        seen = set()
        for l in logs:
            if l in seen: continue
            seen.add(l); print(l[:600])


asyncio.run(main())
