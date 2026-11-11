const assert = require('node:assert/strict');
const {test} = require('node:test');
const {stalled, recover} = require('./runner-fallback.cjs');

const NOW = Date.parse('2026-11-11T16:00:00Z');
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

function fixture({attempt = 1, jobs = [queued], superseded = false, cancelled = true} = {}) {
  const calls = [];
  const run = {id: 123, run_attempt: attempt, status: 'in_progress',
    head_branch: 'topic', event: 'pull_request', created_at: queued.created_at,
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
  return {calls, options: {
    github: {rest: {actions}, async paginate(method, args) {
      if (method === actions.listWorkflowRuns) {
        if (args.workflow_id !== 'check.yml') return [];
        return superseded ? [run, {...run, id: 124, status: 'completed',
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
