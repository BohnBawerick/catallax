// Run from the repository root. Keep Cargo artifacts outside the checkout.
// Usage: CARGO_TARGET_DIR=/tmp/catallax-target nice -n 10 node site/render/collect-evidence.mjs
import { spawnSync } from 'node:child_process';
import { readFileSync, writeFileSync } from 'node:fs';

const root = process.cwd();
const run = (command, args) => {
  console.log(`$ ${command} ${args.join(' ')}`);
  const result = spawnSync(command, args, { encoding: 'utf8', maxBuffer: 64 * 1024 * 1024, env: { ...process.env, CARGO_TERM_COLOR: 'never' } });
  if (result.error || result.status !== 0) throw result.error || new Error(result.stderr + result.stdout);
  return result;
};
const clean = text => text.replaceAll(root, '.').replaceAll(process.env.CARGO_TARGET_DIR || `${root}/target`, 'target').replace(/\x1b\[[0-9;]*m/g, '');
const revision = run('git', ['rev-parse', 'HEAD']).stdout.trim();
const files = run('git', ['ls-files', '-z', '--', 'crates']).stdout.split('\0').filter(Boolean);
const lines = text => text ? text.split('\n').length - (text.endsWith('\n') ? 1 : 0) : 0;
const rust = files.filter(file => file.endsWith('.rs'));
const cpp = files.filter(file => /\.(?:c|cc|cpp|cxx|h|hh|hpp|hxx)$/i.test(file));
const rustText = rust.map(file => readFileSync(file, 'utf8'));
const counts = {
  rust_lines: rustText.reduce((n, text) => n + lines(text), 0),
  cpp_lines: cpp.reduce((n, file) => n + lines(readFileSync(file, 'utf8')), 0),
  unsafe_blocks: rustText.reduce((n, text) => n + [...text.matchAll(/\bunsafe\s*\{/g)].length, 0),
};
const commands = [
  ['cargo', ['build', '--workspace', '--release', '--locked', '-j', '2']],
  ['cargo', ['test', '-p', 'catallax-numbers', '--release', '--lib', '--locked', '-j', '2', '--', '--test-threads=2']],
];
const transcript = commands.map(([command, args]) => {
  const result = run(command, args);
  const output = clean(result.stderr + result.stdout);
  console.log(output);
  return { command: `${command} ${args.join(' ')}`, exit_code: result.status, output };
});
writeFileSync('site/render/terminal-output.json', JSON.stringify({ revision, note: 'Actual output. Paths made repository-relative. Video excerpts are edited for time; this is not a full-suite run.', commands: transcript }, null, 2) + '\n');
const listing = run('cargo', ['nextest', 'list', '--workspace', '--release', '--locked', '--build-jobs', '2', '--message-format', 'json']);
const suites = Object.values(JSON.parse(listing.stdout)['rust-suites']);
const tests = suites.flatMap(suite => Object.values(suite.testcases));
counts.tests = tests.filter(test => !test.ignored).length;
counts.ignored_tests = tests.filter(test => test.ignored).length;
const evidence = {
  revision, counts,
  method: {
    files: 'git ls-files -z -- crates. Only tracked workspace sources; dependencies and generated artifacts excluded.',
    rust_lines: 'Physical lines in .rs files, including comments, blanks and tests. Final newline does not add an extra line.',
    cpp_lines: 'Physical lines in .c, .cc, .cpp, .cxx, .h, .hh, .hpp and .hxx files.',
    unsafe_blocks: 'Matches of the Rust unsafe keyword followed by an opening brace in those .rs files. Zero matches; the workspace forbids unsafe code.',
    tests: 'cargo nextest list --workspace --release --locked --build-jobs 2 --message-format json; testcases with ignored=false. Includes each compiled binary target; excludes doctests.',
  },
  sources: {
    current_implementation: 'https://github.com/BohnBawerick/catallax/blob/' + revision + '/README.md',
    direction: 'https://github.com/BohnBawerick/catallax/blob/' + revision + '/docs/PRD.md',
    ci: 'https://github.com/BohnBawerick/catallax/blob/' + revision + '/.github/workflows/check.yml',
    bio: 'https://pedro.zigrawbohnbawerick.org/',
  },
  historical_robustness: {
    source: 'Project-owner supplied report of the latest kernel, matched against our own hard-test corpus.',
    scope: 'Historical source-run figures, not a measurement of this public revision. The measured source changes are not in the public starting tree. Raw source-run artifacts are not included here; these totals are owner-reported, not independently reproduced by this site capture.',
    requests: { total: 379, answered: 363, refused: 16, impossible_geometry: 6, timeouts: 10 },
    failure_packs: { total: 98, answered: 94 },
    previously_answered_output_changes: 0,
    speed_comparisons: 'Not measured. No cross-kernel speed comparison is claimed.',
  },
};
writeFileSync('site/evidence.json', JSON.stringify(evidence, null, 2) + '\n');
let html = readFileSync('site/index.html', 'utf8');
for (const [key, value] of Object.entries(counts)) {
  html = html.replace(new RegExp(`(data-metric="${key}">)[^<]*`), (_, prefix) => `${prefix}${value.toLocaleString('en-US')}`);
}
const committed = new Date(run('git', ['show', '-s', '--format=%cI', revision]).stdout.trim().slice(0, 10));
const day = new Intl.DateTimeFormat('en-GB', { dateStyle: 'long', timeZone: 'UTC' }).format(committed);
html = html.replace(/(<span data-snapshot>)[^<]*/, `$1${revision.slice(0, 7)}, ${day}`);
writeFileSync('site/index.html', html);
console.log(JSON.stringify(evidence, null, 2));
