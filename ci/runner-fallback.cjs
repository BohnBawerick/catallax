// Only Actions metadata is read here. Never execute a candidate's code with this token.
const OPTIONAL = new Set(['self-hosted', 'boat', 'boat-large']);

function superseded(runs, run) {
  return runs.some(other => other.head_branch === run.head_branch && other.event === run.event
    && (Date.parse(other.created_at) > Date.parse(run.created_at) || other.id > run.id));
}

function stalled(job, now) {
  const optional = (job.labels || []).some(label => OPTIONAL.has(label.toLowerCase()));
  if (!optional || job.runner_id) return false;
  if (job.status === 'completed') return job.conclusion === 'failure';
  const created = Date.parse(job.created_at);
  return job.status === 'queued' && Number.isFinite(created) && now - created >= 300_000;
}

async function recover({github, context, core, now = Date.now, sleep = ms => new Promise(r => setTimeout(r, ms))}) {
  const repo = context.repo;
  // Scheduled sweeps catch late jobs too, including aggregation after a long replay.
  for (const workflow_id of ['check.yml', 'rounded-step-speed.yml']) {
    const runs = await github.paginate(github.rest.actions.listWorkflowRuns, {
      ...repo, workflow_id, per_page: 100,
      created: `>=${new Date(now() - 24 * 60 * 60 * 1000).toISOString()}`,
    });
    for (const run of runs) {
      if (run.run_attempt !== 1 || !['queued', 'in_progress'].includes(run.status)) continue;
      const jobs = await github.paginate(github.rest.actions.listJobsForWorkflowRunAttempt, {
        ...repo, run_id: run.id, attempt_number: run.run_attempt, per_page: 100,
      });
      const stuck = jobs.filter(job => stalled(job, now()));
      if (!stuck.length) continue;
      // A newer push supersedes this run. Do not undo ordinary concurrency cancellation.
      if (superseded(runs, run)) continue;
      const current = (await github.rest.actions.getWorkflowRun({...repo, run_id: run.id})).data;
      if (current.run_attempt !== 1 || !['queued', 'in_progress'].includes(current.status)) continue;
      core.warning(`Runner infrastructure stalled in ${run.html_url}; restarting the whole run on GitHub.`);
      await github.rest.actions.cancelWorkflowRun({...repo, run_id: run.id});
      let cancelled = false;
      for (let attempt = 0; attempt < 12; attempt++) {
        await sleep(5000);
        const state = (await github.rest.actions.getWorkflowRun({...repo, run_id: run.id})).data;
        if (state.run_attempt !== 1) break;
        if (state.status === 'completed') {
          cancelled = state.conclusion === 'cancelled';
          break;
        }
        // An unresponsive runner can delay ordinary cancellation indefinitely.
        if (attempt === 5) await github.rest.actions.forceCancelWorkflowRun({...repo, run_id: run.id});
      }
      if (!cancelled) {
        core.warning(`Cancellation has not completed for ${run.html_url}; no partial rerun was requested.`);
        continue;
      }
      // Cancellation can take a minute. Recheck both runs and the live head;
      // a new push may not have created its workflow run yet.
      const refreshed = await github.paginate(github.rest.actions.listWorkflowRuns, {
        ...repo, workflow_id, branch: run.head_branch, event: run.event, per_page: 100,
        created: `>=${run.created_at}`,
      });
      if (superseded(refreshed, run)) continue;
      let head;
      if (run.event === 'pull_request') {
        const pull_number = run.pull_requests?.[0]?.number;
        if (!pull_number) {
          core.warning(`Cannot confirm the current pull request head for ${run.html_url}; refusing recovery.`);
          continue;
        }
        const pull = (await github.rest.pulls.get({...repo, pull_number})).data;
        if (pull.state !== 'open') continue;
        head = pull.head.sha;
      } else {
        head = (await github.rest.git.getRef({...repo, ref: `heads/${run.head_branch}`})).data.object.sha;
      }
      if (head !== run.head_sha) continue;
      const final = (await github.rest.actions.getWorkflowRun({...repo, run_id: run.id})).data;
      if (final.run_attempt !== 1 || final.status !== 'completed' || final.conclusion !== 'cancelled') continue;
      await github.rest.actions.reRunWorkflow({...repo, run_id: run.id});
      await core.summary.addRaw(`Runner infrastructure fallback: full rerun of ${run.html_url} on GitHub.\n`).write();
    }
  }
}

module.exports = {stalled, recover};
