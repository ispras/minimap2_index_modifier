"""Collinear PAF selection."""
import sys, re
CS = re.compile('(:[0-9]+|\\*[a-z][a-z]|[=\\+\\-][A-Za-z]+|~[a-z]{2}[0-9]+[a-z]{2})')
MAXTRIM = 0.05

def _toklen(tok):
    c = tok[0]
    if c == ':':
        L = int(tok[1:])
        return (L, L)
    if c == '*':
        return (1, 1)
    if c == '=':
        L = len(tok) - 1
        return (L, L)
    if c == '+':
        return (len(tok) - 1, 0)
    if c == '-':
        return (0, len(tok) - 1)
    sys.exit('unexpected cs token %s' % tok)

def _trim_left(toks, q, t, qmin, tmin):
    i = 0
    while i < len(toks) and (q < qmin or t < tmin):
        dq, dt = _toklen(toks[i])
        if toks[i][0] in ':=' and q + dq > qmin and (t + dt > tmin):
            need = max(qmin - q, tmin - t)
            toks[i] = ':%d' % (dq - need)
            q += need
            t += need
            break
        q += dq
        t += dt
        i += 1
    return (toks[i:], q, t)

def trim(r, qs, ts, qe, te):
    toks = CS.findall(r['cs'])
    toks, q, t = _trim_left(toks, r['qs'], r['ts'], qs, ts)
    rtoks = toks[::-1]
    rtoks, rq, rt = _trim_left(rtoks, -r['qe'], -r['te'], -qe, -te)
    toks = rtoks[::-1]
    if not toks or toks[0][0] not in ':=' or toks[-1][0] not in ':=':
        while toks and toks[0][0] not in ':=':
            dq, dt = _toklen(toks[0])
            q += dq
            t += dt
            toks = toks[1:]
        while toks and toks[-1][0] not in ':=':
            dq, dt = _toklen(toks[-1])
            rq += dq
            rt += dt
            toks = toks[:-1]
    if not toks:
        return False
    r['qs'], r['ts'], r['qe'], r['te'], r['cs'] = (q, t, -rq, -rt, ''.join(toks))
    return True

def read_paf(path):
    recs = []
    for n, line in enumerate(open(path), 1):
        f = line.rstrip('\n').split('\t')
        tags = dict(((x[:2], x[5:]) for x in f[12:]))
        recs.append(dict(qn=f[0], ql=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4], tn=f[5], tl=int(f[6]), ts=int(f[7]), te=int(f[8]), nm=int(f[9]), tp=tags.get('tp'), cs=tags.get('cs'), line=n))
    return recs

def select(recs):
    st = dict(total=len(recs), tp_other=0, rev=0, over_q=0, over_t=0, noncol=0, trimmed=0, trim_bp=0)
    cand = []
    for r in recs:
        if r['tp'] not in ('P', 'S'):
            st['tp_other'] += 1
            continue
        if r['strand'] != '+':
            st['rev'] += 1
            continue
        if r['cs'] is None:
            sys.exit('no cs tag on line %d' % r['line'])
        cand.append(r)
    cand.sort(key=lambda r: (-r['nm'], r['qs']))
    st['trimmed'] = 0
    st['trim_bp'] = 0
    kept = []
    for r in cand:
        lim = MAXTRIM * (r['qe'] - r['qs'])
        qs, qe, ts, te = (r['qs'], r['qe'], r['ts'], r['te'])
        bad = False
        for k in kept:
            oq = min(qe, k['qe']) - max(qs, k['qs'])
            ot = min(te, k['te']) - max(ts, k['ts'])
            if oq > lim:
                st['over_q'] += 1
                bad = True
                break
            if ot > lim:
                st['over_t'] += 1
                bad = True
                break
            before = k['qs'] + k['qe'] < qs + qe
            if before != (k['ts'] + k['te'] < ts + te):
                st['noncol'] += 1
                bad = True
                break
            if before:
                qs = max(qs, k['qe'])
                ts = max(ts, k['te'])
            else:
                qe = min(qe, k['qs'])
                te = min(te, k['ts'])
        if bad:
            continue
        if qs >= qe or ts >= te:
            st['over_q'] += 1
            continue
        if (qs, qe, ts, te) != (r['qs'], r['qe'], r['ts'], r['te']):
            n0 = r['qe'] - r['qs']
            if not trim(r, qs, ts, qe, te):
                st['over_q'] += 1
                continue
            st['trimmed'] += 1
            st['trim_bp'] += n0 - (r['qe'] - r['qs'])
        if any((not (r['qe'] <= k['qs'] or r['qs'] >= k['qe']) for k in kept)):
            st['over_q'] += 1
            continue
        if any((not (r['te'] <= k['ts'] or r['ts'] >= k['te']) for k in kept)):
            st['over_t'] += 1
            continue
        if any(((r['qs'] >= k['qe']) != (r['ts'] >= k['te']) for k in kept)):
            st['noncol'] += 1
            continue
        kept.append(r)
    kept.sort(key=lambda r: r['qs'])
    st['kept'] = len(kept)
    return (kept, st)

def fmt(st):
    return 'PAF records %(total)d; not P/S %(tp_other)d; reverse strand %(rev)d; query overlap %(over_q)d, target overlap %(over_t)d; non-collinear %(noncol)d; trimmed %(trimmed)d (%(trim_bp)d bp); kept %(kept)d' % st
