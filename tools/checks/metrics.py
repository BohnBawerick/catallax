#!/usr/bin/env python3
"""Record public check timings and selection decisions after the final gate."""
import argparse
import base64
from datetime import datetime
import json
import os
from pathlib import Path
import subprocess

REPOSITORY = 'BohnBawerick/catallax'
GITHUB_CLI = os.environ.get('CATALLAX_GITHUB_CLI', 'gh-axi')


def api(path):
    # gh-axi's jq string output avoids its human-readable table encoding.
    command = [GITHUB_CLI, 'api', path, '--jq', 'tojson | @base64']
    if GITHUB_CLI == 'gh-axi':
        command.append('--full')
    encoded = subprocess.check_output(command, text=True)
    return json.loads(base64.b64decode(encoded))


def milliseconds(start, end):
    return round((datetime.fromisoformat(end.replace('Z', '+00:00')) -
                  datetime.fromisoformat(start.replace('Z', '+00:00'))).total_seconds() * 1000)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-id', required=True)
    parser.add_argument('--label', required=True)
    parser.add_argument('--selection', required=True)
    parser.add_argument('--replay-summary')
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    run = api(f'repos/{REPOSITORY}/actions/runs/{args.run_id}')
    jobs, page = [], 1
    while True:
        batch = api(f'repos/{REPOSITORY}/actions/runs/{args.run_id}/attempts/{run["run_attempt"]}/jobs?per_page=100&page={page}')['jobs']
        jobs.extend(batch)
        if len(batch) < 100:
            break
        page += 1
    selection = json.loads(Path(args.selection).read_text())
    replay = json.loads(Path(args.replay_summary).read_text()) if args.replay_summary else {
        'selected_cases': 0, 'elapsed_ms': 0}
    if replay['selected_cases'] != selection['selection']['selected_cases']:
        raise ValueError('missing or wrong replay timing aggregate')
    stale = replay.get('commit', selection['identity']['commit']) != selection['identity']['commit']
    records = []
    for job in jobs:
        if job['name'] == 'timings':
            continue  # This observational job cannot measure its own completion.
        queued, start, end = job.get('created_at'), job['started_at'], job['completed_at']
        records.append({'job': job['name'], 'conclusion': job['conclusion'],
                        'queue_ms': milliseconds(queued, start) if queued and start else None,
                        'queue_basis': 'job.created_at' if queued else 'unavailable in GitHub API',
                        'wall_ms': milliseconds(start, end) if start and end else 0})
    gates = [job for job in jobs if job['name'] == 'check']
    if len(gates) != 1 or not gates[0]['completed_at']:
        raise ValueError('timings require exactly one completed final gate')
    result = {'format': 'catallax-check-timing/1', 'label': args.label,
              'run_url': run['html_url'], 'commit': selection['identity']['commit'],
              'head': selection['head'], 'run_attempt': run['run_attempt'],
              'rerun_count': run['run_attempt'] - 1, 'jobs': records,
              'time_to_gate_ms': milliseconds(run['created_at'], gates[0]['completed_at']),
              'selected_cases': replay['selected_cases'],
              'replay_elapsed_ms': 0 if stale else replay['elapsed_ms'],
              'replay_kind': 'stale-restored' if stale else 'fresh' if replay['selected_cases'] else 'not-selected',
              'restored_replay_elapsed_ms': replay['elapsed_ms'] if stale else None,
              'decisions': selection['decisions']}
    Path(args.out).write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(f'{args.label}: {run["html_url"]}; {len(jobs)} jobs; {result["selected_cases"]} cases')


if __name__ == '__main__':
    main()
