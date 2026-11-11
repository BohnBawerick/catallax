#!/usr/bin/env python3
"""Choose optional runners only while their health leases are fresh."""
import json
import os
from pathlib import Path
import time


def fresh(value, now):
    try:
        return 0 <= now - int(value) <= 120
    except (ValueError, TypeError):
        return False


def route(env, now):
    trusted = env.get('TRUSTED_EVENT') == 'true'
    first_attempt = env.get('GITHUB_RUN_ATTEMPT', '1') == '1'
    local = (trusted and first_attempt and env.get('USE_LOCAL_RUNNERS') == 'true'
             and fresh(env.get('LOCAL_OK_AT'), now))
    boat = (trusted and first_attempt and env.get('BOAT_ALLOWED') == 'true'
            and env.get('USE_BOAT') == 'true' and fresh(env.get('BOAT_OK_AT'), now))
    return {'light': ['self-hosted', 'linux', 'x64', 'ci-linux-local'] if local
            else ['ubuntu-latest'], 'boat': str(boat).lower()}


if __name__ == '__main__':
    decision = route(os.environ, time.time())
    print(json.dumps(decision, sort_keys=True))
    with Path(os.environ['GITHUB_OUTPUT']).open('a') as output:
        output.write(f'light={json.dumps(decision["light"])}\nboat={decision["boat"]}\n')
