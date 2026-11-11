const assert = require('node:assert/strict');
const {test} = require('node:test');
const {stalled, recover} = require('./runner-fallback.cjs');

const NOW = Date.parse('2026-11-11T16:00:00Z');
const REVISION = 'a'.repeat(40);
const LEGACY_REVISION = 'b'.repeat(40);
const CAPABILITY = '# catallax-runner-recovery: v1';
const queued = {status: 'queued', labels: ['self-hosted', 'ci-linux-local'],
  runner_id: 0, created_at: '2026-11-11T15:50:00Z'};

test('only optional, unassigned jobs beyond the queue grace period stall', () => {
  assert.equal(stalled(queued, NOW), true);
  assert.equal(stalled({...queued, labels: ['boat-large']}, NOW), true);
  for (const changes of [{labels: ['ubuntu-latest']}, {runner_id: 42},
    {status: 'in_progress'}, {created_at: '2026-11-11T15:59:00Z'}, {created_at: null},
    {status: 'completed', conclusion: 'cancelled'}]) {
    assert.equal(stalled({...queued, ...changes}, NOW), false);
  }
  assert.equal(stalled({...queued, status: 'completed', conclusion: 'failure'}, NOW), true);
});

function fixture({attempt = 1, jobs = [queued], superseded = false,
  supersededDuringCancel = false, headChangesDuringCancel = false,
  event = 'pull_request', cancelled = true, recoveryCapable = true} = {}) {
  const calls = [];
  const run = {id: 123, run_attempt: attempt, status: 'in_progress',
    head_branch: 'topic', head_sha: 'candidate', pull_requests: [{number: 1}],
    event, created_at: queued.created_at, referenced_workflows: [{
      path: `example/project/.github/workflows/runner-route.yml@${REVISION}`,
      ref: event === 'pull_request' ? 'refs/pull/1/merge' : 'refs/heads/topic', sha: REVISION,
    }],
    html_url: 'https://github.com/example/project/actions/runs/123'};
  const actions = {
    listWorkflowRuns() {}, listJobsForWorkflowRunAttempt() {},
    async getWorkflowRun() {
      return {data: calls.includes('cancel') && cancelled
        ? {...run, status: 'completed', conclusion: 'cancelled'} : run};
    },
    async cancelWorkflowRun() { calls.push('cancel'); },
    async forceCancelWorkflowRun() { calls.push('force-cancel'); },
    async reRunWorkflow() { calls.push('rerun-full'); },
  };
  return {calls, run, options: {
    github: {rest: {actions, repos: {async getContent({ref}) {
      return {data: {type: 'file', encoding: 'base64', content: Buffer.from(
        recoveryCapable && ref === REVISION ? CAPABILITY + '\nname: check\n' : 'name: check\n'
      ).toString('base64')}};
    }}, pulls: {async get() {
      return {data: {state: 'open', head: {sha: headChangesDuringCancel ? 'new-push' : run.head_sha}}};
    }}, git: {async getRef() {
      return {data: {object: {sha: headChangesDuringCancel ? 'new-push' : run.head_sha}}};
    }}}, async paginate(method, args) {
      if (method === actions.listWorkflowRuns) {
        if (args.workflow_id !== 'check.yml') return [];
        return (superseded || (supersededDuringCancel && calls.includes('cancel')))
          ? [run, {...run, id: 124, status: 'completed',
          conclusion: 'success', created_at: '2026-11-11T15:59:00Z'}] : [run];
      }
      return jobs;
    }},
    context: {repo: {owner: 'example', repo: 'project'}},
    core: {warning() {}, summary: {addRaw() {return this;}, async write() {}}},
    now: () => NOW, sleep: async () => {},
  }};
}

test('stalled runner cancels and reruns the whole workflow once', async () => {
  const {calls, options} = fixture();
  await recover(options);
  assert.deepEqual(calls, ['cancel', 'rerun-full']);
});

test('reruns, superseded runs and healthy jobs are left alone', async () => {
  for (const input of [{attempt: 2}, {superseded: true}, {jobs: []},
    {jobs: [{...queued, labels: ['ubuntu-latest']}]}, {jobs: [{...queued, runner_id: 42}]}]) {
    const {calls, options} = fixture(input);
    await recover(options);
    assert.deepEqual(calls, []);
  }
});

test('never rerun before cancellation completes', async () => {
  const {calls, options} = fixture({cancelled: false});
  await recover(options);
  assert.deepEqual(calls, ['cancel', 'force-cancel']);
});

test('a newer push during cancellation prevents an obsolete restart', async () => {
  const {calls, options} = fixture({supersededDuringCancel: true});
  await recover(options);
  assert.deepEqual(calls, ['cancel']);
});

test('a new head without a workflow run also prevents an obsolete restart', async () => {
  const {calls, options} = fixture({headChangesDuringCancel: true});
  await recover(options);
  assert.deepEqual(calls, ['cancel']);
});

test('push recovery checks the live branch head as well', async () => {
  for (const changed of [false, true]) {
    const {calls, options} = fixture({event: 'push', headChangesDuringCancel: changed});
    await recover(options);
    assert.deepEqual(calls, changed ? ['cancel'] : ['cancel', 'rerun-full']);
  }
});

test('legacy workflow revisions are left untouched', async () => {
  const {calls, options} = fixture({recoveryCapable: false});
  await recover(options);
  assert.deepEqual(calls, []);
});

test('capability belongs to the recorded workflow revision, not the candidate head', async () => {
  const {calls, run, options} = fixture();
  run.head_sha = REVISION;
  run.referenced_workflows[0].sha = LEGACY_REVISION;
  run.referenced_workflows[0].path = `example/project/.github/workflows/runner-route.yml@${LEGACY_REVISION}`;
  await recover(options);
  assert.deepEqual(calls, []);
});

test('missing revision metadata and unreadable capability fail closed before cancellation', async () => {
  for (const fault of ['missing', 'different-ref', 'unreadable']) {
    const {calls, run, options} = fixture();
    if (fault === 'missing') run.referenced_workflows = [];
    if (fault === 'different-ref') run.referenced_workflows[0].ref = 'refs/heads/main';
    if (fault === 'unreadable') options.github.rest.repos.getContent = async () => {throw new Error('unavailable');};
    await recover(options);
    assert.deepEqual(calls, []);
  }
});

test('a mixed-version sweep recovers only the capable revision', async () => {
  const {calls, run, options} = fixture();
  const legacy = {...run, id: 122, head_branch: 'legacy', pull_requests: [{number: 2}],
    referenced_workflows: [{path: `example/project/.github/workflows/runner-route.yml@${LEGACY_REVISION}`,
      ref: 'refs/pull/2/merge', sha: LEGACY_REVISION}]};
  const {actions} = options.github.rest;
  const paginate = options.github.paginate;
  options.github.paginate = async (method, args) => method === actions.listWorkflowRuns
    ? (args.workflow_id === 'check.yml' ? [legacy, run] : []) : paginate(method, args);
  const cancelled = new Set();
  actions.getWorkflowRun = async ({run_id}) => ({data: {
    ...(run_id === legacy.id ? legacy : run),
    ...(cancelled.has(run_id) ? {status: 'completed', conclusion: 'cancelled'} : {}),
  }});
  actions.cancelWorkflowRun = async ({run_id}) => {cancelled.add(run_id); calls.push(`cancel:${run_id}`);};
  actions.reRunWorkflow = async ({run_id}) => {calls.push(`rerun-full:${run_id}`);};
  await recover(options);
  assert.deepEqual(calls, ['cancel:123', 'rerun-full:123']);
});
