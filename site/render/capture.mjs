// Frame-exact canvas capture through chrome-devtools-axi; ffmpeg encodes H.264.
// Serve site/ first. Usage: SITE_URL=http://127.0.0.1:8873 node site/render/capture.mjs video-a.html site/media/kernel-assembly.mp4
// Use a task-owned CHROME_DEVTOOLS_AXI_SESSION. No other browser is closed.
import { spawn, spawnSync } from 'node:child_process';

const [page, out, crf = '24'] = process.argv.slice(2);
if (!page || !out || !/^[a-z0-9-]+\.html$/.test(page)) throw new Error('Expected render page and output path');
const base = process.env.SITE_URL || 'http://127.0.0.1:8873';
const run = script => {
  const result = spawnSync('chrome-devtools-axi', ['run'], { input: script, encoding: 'utf8', maxBuffer: 64 * 1024 * 1024 });
  if (result.error || result.status !== 0) throw result.error || new Error(result.stderr + result.stdout);
  return result.stdout.trim();
};
const init = run(`await page.open(${JSON.stringify(`${base}/render/${page}`)}); console.log(await page.eval(async () => { await window.__ready; return window.__frames; }));`);
const frames = Number(init);
if (!Number.isSafeInteger(frames) || frames < 1) throw new Error(`Invalid frame count: ${init}`);
const ff = spawn('ffmpeg', ['-y', '-loglevel', 'error', '-f', 'image2pipe', '-framerate', '30', '-c:v', 'png', '-i', '-', '-c:v', 'libx264', '-threads', '2', '-preset', 'slow', '-tune', 'animation', '-crf', crf, '-pix_fmt', 'yuv420p', '-movflags', '+faststart', '-an', out], { stdio: ['pipe', 'inherit', 'inherit'] });
const finished = new Promise((resolve, reject) => {
  ff.once('error', reject);
  ff.once('close', code => code === 0 ? resolve() : reject(new Error(`ffmpeg exit ${code}`)));
});
try {
  for (let start = 0; start < frames; start += 30) {
    const count = Math.min(30, frames - start);
    const data = run(`for (let i = ${start}; i < ${start + count}; i++) console.log(await page.eval('(() => { renderFrame(' + i + '); return document.querySelector("canvas").toDataURL("image/png"); })()'));`);
    const images = data.split('\n').filter(Boolean);
    if (images.length !== count || images.some(image => !image.startsWith('data:image/png;base64,'))) throw new Error(`Missing frames at ${start}`);
    for (const image of images) {
      if (!ff.stdin.write(Buffer.from(image.slice('data:image/png;base64,'.length), 'base64'))) await new Promise(resolve => ff.stdin.once('drain', resolve));
    }
  }
  ff.stdin.end();
  await finished;
} catch (error) {
  ff.kill();
  await finished.catch(() => {});
  throw error;
}
console.log(`${out}: ${frames} frames`);
