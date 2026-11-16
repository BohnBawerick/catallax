export async function measure(wasm, body, iterations = 7) {
  const { instance } = await WebAssembly.instantiate(wasm, {});
  const e = instance.exports;
  const buffer = e.input_buffer(body.byteLength);
  new Uint8Array(e.memory.buffer, buffer, body.byteLength).set(new Uint8Array(body));
  let start = performance.now();
  const cold_read_work = e.probe_read();
  const cold_read_ms = performance.now() - start;
  e.source_from_input();
  e.clear();
  for (let i = 0; i < 3; i++) {
    e.probe_read();
    e.clear();
  }
  const rows = [];
  for (let i = 0; i < iterations; i++) {
    e.clear();
    start = performance.now();
    const write_work = e.probe_write();
    const write_ms = performance.now() - start;
    start = performance.now();
    const read_work = e.probe_read();
    const read_ms = performance.now() - start;
    rows.push({ write_ms, read_ms, write_work, read_work });
    await new Promise(resolve => setTimeout(resolve, 0));
  }
  return { bytes: body.byteLength, cold_read_ms, cold_read_work, rows };
}

if (typeof process !== 'undefined' && process.versions?.node) {
  const { readFile } = await import('node:fs/promises');
  const [wasm, body] = await Promise.all(process.argv.slice(2, 4).map(path => readFile(path)));
  console.log(JSON.stringify(await measure(wasm, body)));
}
