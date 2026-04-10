const IUPAC = {
  A: 'A', C: 'C', G: 'G', T: 'T', U: 'T',
  R: 'AG', Y: 'CT', S: 'GC', W: 'AT', K: 'GT', M: 'AC',
  B: 'CGT', D: 'AGT', H: 'ACT', V: 'ACG', N: 'ACGT'
};

const state = {
  analyses: { seq1: null, seq2: null },
  fragmentLookup: new Map(),
  syntheticBlocks: [],
  assembly: []
};

function $(id) { return document.getElementById(id); }
function uid(prefix = 'id') { return `${prefix}-${Math.random().toString(36).slice(2, 10)}`; }

function revComp(seq) {
  const map = { A: 'T', T: 'A', G: 'C', C: 'G', R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K', B: 'V', D: 'H', H: 'D', V: 'B', N: 'N' };
  return seq.toUpperCase().split('').reverse().map(ch => map[ch] || 'N').join('');
}

function baseMatches(base, code) {
  return (IUPAC[code] || 'ACGT').includes(base);
}

function cleanSequence(raw) {
  if (!raw) return '';
  const text = raw.trim();
  if (!text) return '';
  if (/^LOCUS\b/m.test(text) && /^ORIGIN\b/m.test(text)) {
    const origin = text.split(/^ORIGIN\b/m)[1] || '';
    const seqChunk = origin.split(/^\/\//m)[0] || origin;
    return (seqChunk.match(/[A-Za-z]/g) || []).join('').toUpperCase().replace(/U/g, 'T').replace(/[^ACGTRYSWKMBDHVN]/g, '');
  }
  return text.split(/\r?\n/).filter(line => !line.startsWith('>')).join('').toUpperCase().replace(/U/g, 'T').replace(/[^ACGTRYSWKMBDHVN]/g, '');
}

function formatCut(enzyme) {
  const c = enzyme.cut;
  if (c.kind === 'simple') return `${enzyme.recognition}(${c.top}/${c.bottom})`;
  return `(${c.leftTop}/${c.leftBottom})${enzyme.recognition}(${c.rightTop}/${c.rightBottom})`;
}

function populateEnzymeSelect(select) {
  select.innerHTML = '';
  window.TYPE_IIS_ENZYMES.slice().sort((a, b) => a.name.localeCompare(b.name)).forEach(enzyme => {
    const opt = document.createElement('option');
    const aliasText = enzyme.aliases?.length ? ` | aliases: ${enzyme.aliases.join(', ')}` : '';
    opt.value = enzyme.name;
    opt.textContent = `${enzyme.name} — ${formatCut(enzyme)}${aliasText}`;
    select.appendChild(opt);
  });
}

function getSelectedEnzymes(selectId) {
  const names = Array.from($(selectId).selectedOptions).map(opt => opt.value);
  return window.TYPE_IIS_ENZYMES.filter(e => names.includes(e.name));
}

function findSites(seq, enzymes) {
  const hits = [];
  const upper = seq.toUpperCase();
  for (const enzyme of enzymes) {
    const site = enzyme.recognition.toUpperCase();
    const L = site.length;
    const rc = revComp(site);
    for (let i = 0; i <= upper.length - L; i += 1) {
      const fragment = upper.slice(i, i + L);
      let forward = true;
      let reverse = true;
      for (let j = 0; j < L; j += 1) {
        if (!baseMatches(fragment[j], site[j])) forward = false;
        if (!baseMatches(fragment[j], rc[j])) reverse = false;
        if (!forward && !reverse) break;
      }
      if (forward) hits.push(buildHit(upper, enzyme, i, 'forward'));
      if (reverse && rc !== site) hits.push(buildHit(upper, enzyme, i, 'reverse'));
    }
  }
  return hits.sort((a, b) => a.start - b.start || a.enzyme.name.localeCompare(b.enzyme.name) || a.orientation.localeCompare(b.orientation));
}

function buildHit(seq, enzyme, start, orientation) {
  const L = enzyme.recognition.length;
  const end = start + L;
  const hit = { id: `${enzyme.name}:${orientation}:${start}`, enzyme, orientation, start, end, recognitionShown: orientation === 'forward' ? enzyme.recognition : revComp(enzyme.recognition), cuts: [] };
  if (enzyme.cut.kind === 'simple') {
    let topCut, bottomCut;
    if (orientation === 'forward') {
      topCut = end + enzyme.cut.top;
      bottomCut = end + enzyme.cut.bottom;
    } else {
      topCut = start - enzyme.cut.bottom;
      bottomCut = start - enzyme.cut.top;
    }
    hit.cuts.push(makeCutDescriptor(seq, topCut, bottomCut, enzyme, orientation, start, 'single'));
  } else {
    let leftTop, leftBottom, rightTop, rightBottom;
    if (orientation === 'forward') {
      leftTop = start - enzyme.cut.leftTop;
      leftBottom = start - enzyme.cut.leftBottom;
      rightTop = end + enzyme.cut.rightTop;
      rightBottom = end + enzyme.cut.rightBottom;
    } else {
      leftTop = start - enzyme.cut.rightBottom;
      leftBottom = start - enzyme.cut.rightTop;
      rightTop = end + enzyme.cut.leftBottom;
      rightBottom = end + enzyme.cut.leftTop;
    }
    hit.cuts.push(makeCutDescriptor(seq, leftTop, leftBottom, enzyme, orientation, start, 'left'));
    hit.cuts.push(makeCutDescriptor(seq, rightTop, rightBottom, enzyme, orientation, start, 'right'));
  }
  return hit;
}

function makeCutDescriptor(seq, topCut, bottomCut, enzyme, orientation, siteStart, flank = 'single') {
  const inBounds = topCut >= 0 && bottomCut >= 0 && topCut <= seq.length && bottomCut <= seq.length;
  const leftEnd = inBounds ? deriveFragmentEnd(seq, topCut, bottomCut, 'left') : null;
  const rightEnd = inBounds ? deriveFragmentEnd(seq, topCut, bottomCut, 'right') : null;
  return { id: `${enzyme.name}:${orientation}:${siteStart}:${flank}:${topCut}:${bottomCut}`, topCut, bottomCut, overhangLength: Math.abs(topCut - bottomCut), flank, inBounds, leftEnd, rightEnd, boundary: Math.min(topCut, bottomCut) };
}

function deriveFragmentEnd(seq, topCut, bottomCut, side) {
  if (topCut === bottomCut) return { type: 'blunt', sequence: '', length: 0, topCut, bottomCut, side };
  if (side === 'left') {
    if (topCut < bottomCut) return { type: '5prime', sequence: seq.slice(topCut, bottomCut), length: bottomCut - topCut, topCut, bottomCut, side };
    return { type: '3prime', sequence: seq.slice(bottomCut, topCut), length: topCut - bottomCut, topCut, bottomCut, side };
  }
  if (topCut < bottomCut) return { type: '3prime', sequence: seq.slice(topCut, bottomCut), length: bottomCut - topCut, topCut, bottomCut, side };
  return { type: '5prime', sequence: seq.slice(bottomCut, topCut), length: topCut - bottomCut, topCut, bottomCut, side };
}

function flipEnd(end) {
  if (!end) return null;
  if (end.type === 'blunt') return { ...end };
  return { ...end, type: end.type === '5prime' ? '3prime' : '5prime', sequence: revComp(end.sequence) };
}

function compareEnds(endA, endB) {
  if (!endA || !endB) return { compatible: false, reason: 'Missing end information.' };
  if (endA.type === 'blunt' && endB.type === 'blunt') return { compatible: true, reason: 'Both ends are blunt.' };
  if (endA.type === 'blunt' || endB.type === 'blunt') return { compatible: false, reason: 'One end is blunt and the other is sticky.' };
  if (endA.type === endB.type) return { compatible: false, reason: 'Sticky ends have the same polarity.' };
  if (endA.length !== endB.length) return { compatible: false, reason: `Overhang lengths differ (${endA.length} vs ${endB.length}).` };
  const expected = revComp(endA.sequence);
  const compatible = expected === endB.sequence;
  return { compatible, reason: compatible ? `Compatible sticky ends: ${endA.sequence} pairs with ${endB.sequence}.` : `Sequences are not complementary: expected ${expected}, found ${endB.sequence}.` };
}

function describeEnd(end) {
  if (!end) return 'n/a';
  if (end.type === 'blunt') return 'blunt';
  return `${end.type === '5prime' ? '5′' : '3′'} ${end.length} nt (${end.sequence})`;
}

function analyzeSequence(which) {
  const seqRaw = $(`${which}Input`).value;
  const name = $(`${which}Name`).value.trim() || which;
  const enzymes = getSelectedEnzymes(which === 'seq1' ? 'enzymes1' : 'enzymes2');
  const seq = cleanSequence(seqRaw);
  if (!seq) { alert(`No valid sequence found in ${name}.`); return null; }
  if (!enzymes.length) { alert('Select at least one Type IIS enzyme first.'); return null; }
  const hits = findSites(seq, enzymes);
  const analysis = { which, name, seq, hits, enzymes, fragments: [] };
  state.analyses[which] = analysis;
  renderAnalysis(which);
  return analysis;
}

function renderAnalysis(which) {
  const analysis = state.analyses[which];
  if (!analysis) return;
  const summary = $(which === 'seq1' ? 'summary1' : 'summary2');
  const viewer = $(which === 'seq1' ? 'viewer1' : 'viewer2');
  const inBoundsCuts = analysis.hits.flatMap(hit => hit.cuts.filter(c => c.inBounds)).length;
  summary.textContent = `${analysis.name}: ${analysis.seq.length} bp, ${analysis.hits.length} recognition-site hit(s), ${inBoundsCuts} in-bounds cut event(s).`;
  if (!analysis.hits.length) {
    viewer.classList.add('empty');
    viewer.textContent = 'No selected Type IIS sites found in this sequence.';
    return;
  }
  viewer.classList.remove('empty');
  viewer.innerHTML = '';
  analysis.hits.forEach(hit => viewer.appendChild(renderHitBlock(analysis.seq, hit)));
}

function renderHitBlock(seq, hit) {
  const wrap = document.createElement('div');
  wrap.className = 'viewer-block';
  const label = document.createElement('div');
  label.className = 'viewer-label';
  label.textContent = `${hit.enzyme.name} | ${hit.orientation} strand match | site ${hit.start + 1}-${hit.end}`;
  wrap.appendChild(label);

  const windowStart = Math.max(0, Math.min(hit.start, ...hit.cuts.map(c => Math.min(c.topCut, c.bottomCut))) - 8);
  const windowEnd = Math.min(seq.length, Math.max(hit.end, ...hit.cuts.map(c => Math.max(c.topCut, c.bottomCut))) + 8);
  const fragment = {
    seq: seq.slice(windowStart, windowEnd),
    leftEnd: { type: 'blunt', sequence: '', length: 0 },
    rightEnd: { type: 'blunt', sequence: '', length: 0 },
    annotations: [{ start: hit.start - windowStart, end: hit.end - windowStart, className: 'site', label: hit.enzyme.name }],
    overhangWindows: hit.cuts.filter(cut => cut.inBounds).map(cut => ({ start: Math.min(cut.topCut, cut.bottomCut) - windowStart, end: Math.max(cut.topCut, cut.bottomCut) - windowStart, className: cut.topCut < cut.bottomCut ? 'overhang5' : 'overhang3' }))
  };
  wrap.appendChild(renderDsFragment(fragment, { truncate: false, showEnds: false }));

  const meta = document.createElement('div');
  meta.className = 'viewer-hit-meta';
  meta.innerHTML = hit.cuts.map(cut => {
    if (!cut.inBounds) return `${cut.flank}: cut outside sequence bounds (top ${cut.topCut}, bottom ${cut.bottomCut})`;
    return `${cut.flank}: top cut ${cut.topCut}, bottom cut ${cut.bottomCut}, left end ${describeEnd(cut.leftEnd)}, right end ${describeEnd(cut.rightEnd)}`;
  }).join('<br>');
  wrap.appendChild(meta);
  return wrap;
}

function collectFragmentAnnotations(analysis, start, end) {
  const annotations = [];
  analysis.hits.forEach(hit => {
    if (hit.start >= start && hit.end <= end) {
      annotations.push({ start: hit.start - start, end: hit.end - start, className: 'site', label: hit.enzyme.name });
    }
  });
  return annotations;
}

function digestSequence(which) {
  const analysis = state.analyses[which] || analyzeSequence(which);
  if (!analysis) return;
  const simpleCuts = analysis.hits.flatMap(hit => hit.cuts.filter(cut => hit.enzyme.cut.kind === 'simple' && cut.inBounds)).sort((a, b) => a.boundary - b.boundary || a.topCut - b.topCut || a.bottomCut - b.bottomCut);
  const uniqueBoundaries = [0];
  simpleCuts.forEach(cut => {
    if (cut.boundary > 0 && cut.boundary < analysis.seq.length && !uniqueBoundaries.includes(cut.boundary)) uniqueBoundaries.push(cut.boundary);
  });
  if (!uniqueBoundaries.includes(analysis.seq.length)) uniqueBoundaries.push(analysis.seq.length);
  uniqueBoundaries.sort((a, b) => a - b);
  const boundaryCuts = new Map(simpleCuts.map(cut => [cut.boundary, cut]));
  analysis.fragments = [];

  for (let i = 0; i < uniqueBoundaries.length - 1; i += 1) {
    const start = uniqueBoundaries[i];
    const end = uniqueBoundaries[i + 1];
    const seq = analysis.seq.slice(start, end);
    if (!seq) continue;
    const fragment = {
      id: uid(`${which}-fragment`),
      sourceKey: which,
      sourceName: analysis.name,
      label: `${analysis.name} F${analysis.fragments.length + 1}`,
      start,
      end,
      seq,
      annotations: collectFragmentAnnotations(analysis, start, end),
      leftEnd: start === 0 ? { type: 'blunt', sequence: '', length: 0, side: 'left' } : (boundaryCuts.get(start)?.rightEnd || { type: 'blunt', sequence: '', length: 0, side: 'left' }),
      rightEnd: end === analysis.seq.length ? { type: 'blunt', sequence: '', length: 0, side: 'right' } : (boundaryCuts.get(end)?.leftEnd || { type: 'blunt', sequence: '', length: 0, side: 'right' })
    };
    state.fragmentLookup.set(fragment.id, fragment);
    analysis.fragments.push(fragment);
  }
  renderFragmentBank(which);
}

function orientationVariant(fragment, orientation = 'forward') {
  if (orientation === 'forward') {
    return { baseId: fragment.id, label: fragment.label, sourceName: fragment.sourceName, orientation: 'forward', seq: fragment.seq, leftEnd: structuredCloneSafe(fragment.leftEnd), rightEnd: structuredCloneSafe(fragment.rightEnd), start: fragment.start ?? 0, end: fragment.end ?? fragment.seq.length, annotations: structuredCloneSafe(fragment.annotations || []) };
  }
  return {
    baseId: fragment.id,
    label: fragment.label,
    sourceName: fragment.sourceName,
    orientation: 'reverse',
    seq: revComp(fragment.seq),
    leftEnd: flipEnd(fragment.rightEnd),
    rightEnd: flipEnd(fragment.leftEnd),
    start: fragment.start ?? 0,
    end: fragment.end ?? fragment.seq.length,
    annotations: (fragment.annotations || []).map(a => ({ ...a, start: fragment.seq.length - a.end, end: fragment.seq.length - a.start }))
  };
}

function structuredCloneSafe(obj) { return obj ? JSON.parse(JSON.stringify(obj)) : obj; }

function renderFragmentBank(which) {
  const bank = $(which === 'seq1' ? 'bank1' : 'bank2');
  const analysis = state.analyses[which];
  bank.innerHTML = '';
  if (!analysis?.fragments?.length) {
    bank.className = 'bank empty-bank';
    bank.textContent = 'No draggable fragments could be generated.';
    return;
  }
  bank.className = 'bank';
  analysis.fragments.forEach(fragment => bank.appendChild(renderFragmentCard(fragment, { mode: 'bank' })));
}

function renderSyntheticBank() {
  const bank = $('builderBank');
  if (!bank) return;
  bank.innerHTML = '';
  if (!state.syntheticBlocks.length) {
    bank.className = 'bank empty-bank';
    bank.textContent = 'Built restriction blocks will appear here as draggable dsDNA parts.';
    return;
  }
  bank.className = 'bank';
  state.syntheticBlocks.forEach(fragment => bank.appendChild(renderFragmentCard(fragment, { mode: 'bank' })));
}

function endClass(end) {
  return end.type === 'blunt' ? 'blunt' : (end.type === '5prime' ? 'end5' : 'end3');
}

function createStickyGraphic(end, side) {
  const wrap = document.createElement('div');
  wrap.className = `sticky-edge ${side} ${endClass(end)}`;
  const notch = document.createElement('div');
  notch.className = 'sticky-notch';
  const text = document.createElement('div');
  text.className = 'sticky-text';
  text.textContent = end.type === 'blunt' ? 'blunt' : `${end.type === '5prime' ? '5′' : '3′'} ${end.sequence}`;
  wrap.appendChild(notch);
  wrap.appendChild(text);
  return wrap;
}

function baseClassAt(index, annotations, overhangWindows) {
  const classes = [];
  annotations.forEach(a => { if (index >= a.start && index < a.end) classes.push(a.className); });
  overhangWindows.forEach(a => { if (index >= a.start && index < a.end) classes.push(a.className); });
  return classes.join(' ');
}

function makeOverhangWindows(fragment) {
  const windows = [];
  if (fragment.leftEnd && fragment.leftEnd.type !== 'blunt' && fragment.leftEnd.length > 0) {
    windows.push({ start: 0, end: Math.min(fragment.leftEnd.length, fragment.seq.length), className: fragment.leftEnd.type === '5prime' ? 'overhang5' : 'overhang3' });
  }
  if (fragment.rightEnd && fragment.rightEnd.type !== 'blunt' && fragment.rightEnd.length > 0) {
    windows.push({ start: Math.max(0, fragment.seq.length - fragment.rightEnd.length), end: fragment.seq.length, className: fragment.rightEnd.type === '5prime' ? 'overhang5' : 'overhang3' });
  }
  return windows;
}

function renderSequenceChunk(fragment, chunkSeq, offset, annotations, overhangWindows) {
  const chunk = document.createElement('div');
  chunk.className = 'dsdna-chunk';
  const top = document.createElement('div');
  top.className = 'strand-line';
  const bottom = document.createElement('div');
  bottom.className = 'strand-line';
  const rc = revComp(chunkSeq);
  for (let i = 0; i < chunkSeq.length; i += 1) {
    const globalIndex = offset + i;
    const klass = baseClassAt(globalIndex, annotations, overhangWindows);
    const topBase = document.createElement('span');
    topBase.className = `base ${klass}`.trim();
    topBase.textContent = chunkSeq[i];
    const bottomBase = document.createElement('span');
    bottomBase.className = `base ${klass}`.trim();
    bottomBase.textContent = rc[chunkSeq.length - 1 - i];
    top.appendChild(topBase);
    bottom.appendChild(bottomBase);
  }
  chunk.appendChild(top);
  chunk.appendChild(bottom);
  return chunk;
}

function renderDsFragment(fragment, options = {}) {
  const container = document.createElement('div');
  container.className = 'dsdna-view';
  const annotations = fragment.annotations || [];
  const overhangWindows = fragment.overhangWindows || makeOverhangWindows(fragment);
  const showEnds = options.showEnds !== false;
  const truncate = options.truncate !== false;
  const displaySeq = truncate && fragment.seq.length > 180 ? `${fragment.seq.slice(0, 180)}` : fragment.seq;
  const chunkSize = options.chunkSize || 60;

  if (showEnds) {
    const endsRow = document.createElement('div');
    endsRow.className = 'ends-row';
    const left = document.createElement('div');
    left.className = `end-tag ${endClass(fragment.leftEnd)}`;
    left.textContent = `Left end: ${describeEnd(fragment.leftEnd)}`;
    const right = document.createElement('div');
    right.className = `end-tag ${endClass(fragment.rightEnd)}`;
    right.textContent = `Right end: ${describeEnd(fragment.rightEnd)}`;
    endsRow.appendChild(left);
    endsRow.appendChild(right);
    container.appendChild(endsRow);
  }

  for (let offset = 0; offset < displaySeq.length; offset += chunkSize) {
    container.appendChild(renderSequenceChunk(fragment, displaySeq.slice(offset, offset + chunkSize), offset, annotations, overhangWindows));
  }
  if (truncate && fragment.seq.length > displaySeq.length) {
    const more = document.createElement('div');
    more.className = 'fragment-meta';
    more.textContent = `Sequence display truncated to ${displaySeq.length} nt for readability. Full fragment length: ${fragment.seq.length} nt.`;
    container.appendChild(more);
  }
  return container;
}

function renderFragmentCard(item, options = {}) {
  const mode = options.mode || 'bank';
  const fragment = mode === 'assembly' ? item.variant : orientationVariant(item, 'forward');
  const card = document.createElement('div');
  card.className = `fragment ${mode === 'assembly' ? 'assembly-fragment' : 'bank-fragment'}`;
  card.draggable = true;
  if (mode === 'bank') {
    card.dataset.dragType = 'bank';
    card.dataset.fragmentId = item.id;
  } else {
    card.dataset.dragType = 'assembly';
    card.dataset.placedId = item.placedId;
  }

  card.appendChild(createStickyGraphic(fragment.leftEnd, 'left'));
  card.appendChild(createStickyGraphic(fragment.rightEnd, 'right'));

  const head = document.createElement('div');
  head.className = 'fragment-header';
  const orientationBadge = mode === 'assembly'
    ? `<span class="orientation-badge ${fragment.orientation === 'reverse' ? 'reverse' : ''}">${fragment.orientation === 'reverse' ? 'auto-flipped' : 'forward'}</span>`
    : '<span class="orientation-badge">part</span>';
  head.innerHTML = `<div><div class="fragment-title">${fragment.label}</div><div class="fragment-meta">${fragment.sourceName} | ${fragment.start}-${fragment.end} | ${fragment.seq.length} bp</div></div>${orientationBadge}`;
  card.appendChild(head);

  card.appendChild(renderDsFragment(fragment, { truncate: mode !== 'assembly', showEnds: true, chunkSize: mode === 'assembly' ? 72 : 48 }));

  if (mode === 'assembly' && item.joinNote) {
    const note = document.createElement('div');
    note.className = `fragment-meta join-note ${item.joinNote.compatible ? 'ok' : 'bad'}`;
    note.textContent = item.joinNote.text;
    card.appendChild(note);
  }

  card.addEventListener('dragstart', ev => {
    card.classList.add('dragging');
    ev.dataTransfer.effectAllowed = 'move';
    if (mode === 'bank') ev.dataTransfer.setData('application/json', JSON.stringify({ type: 'bank', fragmentId: item.id }));
    else ev.dataTransfer.setData('application/json', JSON.stringify({ type: 'assembly', placedId: item.placedId }));
  });
  card.addEventListener('dragend', () => card.classList.remove('dragging'));
  return card;
}

function validatePlacement(baseFragmentId, slotIndex, movingPlacedId = null) {
  const baseFragment = state.fragmentLookup.get(baseFragmentId);
  if (!baseFragment) return { allowed: false, reason: 'Unknown fragment.' };
  const simulated = state.assembly.filter(item => item.placedId !== movingPlacedId);
  const leftNeighbor = simulated[slotIndex - 1] || null;
  const rightNeighbor = simulated[slotIndex] || null;
  const variants = [orientationVariant(baseFragment, 'forward'), orientationVariant(baseFragment, 'reverse')];
  const evaluations = variants.map(variant => {
    let leftOk = true; let rightOk = true; let leftCmp = null; let rightCmp = null;
    if (leftNeighbor) { leftCmp = compareEnds(leftNeighbor.variant.rightEnd, variant.leftEnd); leftOk = leftCmp.compatible; }
    if (rightNeighbor) { rightCmp = compareEnds(variant.rightEnd, rightNeighbor.variant.leftEnd); rightOk = rightCmp.compatible; }
    return { variant, allowed: leftOk && rightOk, score: Number(leftOk) + Number(rightOk) + (variant.orientation === 'reverse' ? 0 : 0.1), leftCmp, rightCmp, reason: [leftCmp && !leftCmp.compatible ? `Left join: ${leftCmp.reason}` : '', rightCmp && !rightCmp.compatible ? `Right join: ${rightCmp.reason}` : ''].filter(Boolean).join(' ') };
  });
  const winner = evaluations.filter(e => e.allowed).sort((a, b) => b.score - a.score)[0];
  if (winner) return { allowed: true, variant: winner.variant, leftCmp: winner.leftCmp, rightCmp: winner.rightCmp };
  return { allowed: false, reason: evaluations[0]?.reason || 'This fragment does not fit this slot.' };
}

function placePayloadInSlot(payload, slotIndex) {
  const currentLength = state.assembly.length;
  const normalizedSlot = Math.max(0, Math.min(slotIndex, currentLength));
  if (payload.type === 'bank') {
    const validation = validatePlacement(payload.fragmentId, normalizedSlot, null);
    if (!validation.allowed) { flashAssemblyMessage(`Rejected placement. ${validation.reason}`); return; }
    state.assembly.splice(normalizedSlot, 0, { placedId: uid('placed'), baseId: payload.fragmentId, variant: validation.variant, joinNote: null });
  } else if (payload.type === 'assembly') {
    const movingIndex = state.assembly.findIndex(item => item.placedId === payload.placedId);
    if (movingIndex < 0) return;
    const moving = state.assembly[movingIndex];
    const baseId = moving.baseId;
    state.assembly.splice(movingIndex, 1);
    const adjustedSlot = movingIndex < normalizedSlot ? normalizedSlot - 1 : normalizedSlot;
    const validation = validatePlacement(baseId, adjustedSlot, payload.placedId);
    if (!validation.allowed) {
      state.assembly.splice(movingIndex, 0, moving);
      flashAssemblyMessage(`Rejected move. ${validation.reason}`);
      renderAssemblyTrack();
      return;
    }
    moving.variant = validation.variant;
    state.assembly.splice(adjustedSlot, 0, moving);
  }
  refreshAssemblyAnnotations();
  renderAssemblyTrack();
}

function refreshAssemblyAnnotations() {
  state.assembly.forEach(item => { item.joinNote = null; });
  for (let i = 1; i < state.assembly.length; i += 1) {
    const left = state.assembly[i - 1];
    const right = state.assembly[i];
    const cmp = compareEnds(left.variant.rightEnd, right.variant.leftEnd);
    right.joinNote = { compatible: cmp.compatible, text: cmp.compatible ? `Snapped to previous fragment. ${cmp.reason}` : `Join failed. ${cmp.reason}` };
  }
}

function renderAssemblyTrack() {
  const track = $('assemblyTrack');
  track.innerHTML = '';
  const slotCount = Math.max(state.assembly.length + 1, 3);
  for (let i = 0; i < slotCount; i += 1) {
    const dropZone = document.createElement('div');
    dropZone.className = 'drop-zone';
    dropZone.dataset.slotIndex = i;
    dropZone.innerHTML = `<div class="drop-zone-label">Snap</div>`;
    attachDropHandlers(dropZone, i);
    track.appendChild(dropZone);
    if (state.assembly[i]) {
      const placedWrap = document.createElement('div');
      placedWrap.className = 'placed-fragment-wrap';
      placedWrap.appendChild(renderFragmentCard(state.assembly[i], { mode: 'assembly' }));
      track.appendChild(placedWrap);
    }
  }
  updateAssemblyStatus();
}

function attachDropHandlers(el, slotIndex) {
  el.addEventListener('dragover', ev => { ev.preventDefault(); el.classList.add('drag-over'); });
  el.addEventListener('dragleave', () => el.classList.remove('drag-over'));
  el.addEventListener('drop', ev => {
    ev.preventDefault();
    el.classList.remove('drag-over');
    try {
      const payload = JSON.parse(ev.dataTransfer.getData('application/json'));
      placePayloadInSlot(payload, slotIndex);
    } catch {
      flashAssemblyMessage('Unable to read the dragged fragment payload.');
    }
  });
}

function updateAssemblyStatus() {
  const status = $('assemblyStatus');
  const assembled = $('assembledSequence');
  if (!state.assembly.length) {
    status.textContent = 'No fragments in the workbench yet.';
    assembled.textContent = 'Assembled sequence will appear here.';
    return;
  }
  const joins = [];
  let allCompatible = true;
  for (let i = 0; i < state.assembly.length - 1; i += 1) {
    const cmp = compareEnds(state.assembly[i].variant.rightEnd, state.assembly[i + 1].variant.leftEnd);
    if (!cmp.compatible) allCompatible = false;
    joins.push(`Join ${i + 1}→${i + 2}: ${cmp.compatible ? 'compatible' : 'not compatible'} (${cmp.reason})`);
  }
  const mergedSeq = state.assembly.map(item => item.variant.seq).join('');
  status.textContent = allCompatible ? `Workbench contains ${state.assembly.length} fragment(s). Every current join is valid and snapped.` : `Workbench contains ${state.assembly.length} fragment(s). At least one join is invalid.`;
  assembled.innerHTML = `<strong>${allCompatible ? 'Assembled sequence' : 'Provisional assembled sequence'}</strong><br>${mergedSeq}<br><br>${joins.join('<br>') || 'Single fragment only.'}`;
}

function flashAssemblyMessage(message) { $('assemblyStatus').textContent = message; }

function buildSyntheticFragment(enzyme, left, payload, right) {
  const recognition = enzyme.recognition.toUpperCase();
  const seq = `${left}${recognition}${payload}${right}`;
  const block = {
    id: uid('synthetic-block'),
    sourceKey: 'builder',
    sourceName: 'Restriction block builder',
    label: `${enzyme.name} block ${state.syntheticBlocks.length + 1}`,
    start: 0,
    end: seq.length,
    seq,
    annotations: [{ start: left.length, end: left.length + recognition.length, className: 'site', label: enzyme.name }],
    leftEnd: { type: 'blunt', sequence: '', length: 0, side: 'left' },
    rightEnd: { type: 'blunt', sequence: '', length: 0, side: 'right' }
  };
  const analysis = { seq, hits: findSites(seq, [enzyme]) };
  const simpleCut = analysis.hits.flatMap(hit => hit.cuts.filter(c => c.inBounds && hit.enzyme.cut.kind === 'simple'))[0];
  if (simpleCut) {
    block.leftEnd = simpleCut.leftEnd;
    block.rightEnd = simpleCut.rightEnd;
  }
  return block;
}

function applyBuilder() {
  const target = $('builderTarget').value;
  const enzyme = window.TYPE_IIS_ENZYMES.find(e => e.name === $('builderEnzyme').value);
  const left = cleanSequence($('builderLeft').value);
  const payload = cleanSequence($('builderPayload').value);
  const right = cleanSequence($('builderRight').value);
  const mode = $('builderMode').value;
  const block = buildSyntheticFragment(enzyme, left, payload, right);
  const current = cleanSequence($(`${target}Input`).value);
  let updated = block.seq;
  if (mode === 'prefix') updated = `${block.seq}${current}`;
  if (mode === 'suffix') updated = `${current}${block.seq}`;
  $(`${target}Input`).value = updated;
  state.syntheticBlocks.unshift(block);
  state.fragmentLookup.set(block.id, block);
  renderSyntheticBank();
  const preview = $('builderPreview');
  preview.innerHTML = '';
  const blurb = document.createElement('div');
  blurb.className = 'builder-preview-text';
  blurb.innerHTML = `Built a real draggable <strong>${enzyme.name}</strong> block and also applied it to ${target === 'seq1' ? 'Sequence 1' : 'Sequence 2'}. Recognition site remains color-coded in the part below.`;
  preview.appendChild(blurb);
  preview.appendChild(renderDsFragment(block, { truncate: false, showEnds: true, chunkSize: 48 }));
}

function copyBlockOnly() {
  const enzyme = window.TYPE_IIS_ENZYMES.find(e => e.name === $('builderEnzyme').value);
  const left = cleanSequence($('builderLeft').value);
  const payload = cleanSequence($('builderPayload').value);
  const right = cleanSequence($('builderRight').value);
  const block = `${left}${enzyme.recognition}${payload}${right}`;
  navigator.clipboard.writeText(block);
  $('builderPreview').innerHTML = `Copied block sequence for <strong>${enzyme.name}</strong>: ${block}`;
}

function readFileIntoInput(fileInputId, targetTextareaId) {
  $(fileInputId).addEventListener('change', ev => {
    const file = ev.target.files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = e => { $(targetTextareaId).value = e.target.result; };
    reader.readAsText(file);
  });
}

function clearSequence(which) {
  $(`${which}Input`).value = '';
  const viewer = $(which === 'seq1' ? 'viewer1' : 'viewer2');
  viewer.textContent = `Analyze ${which === 'seq1' ? 'sequence 1' : 'sequence 2'} to see the visual digest.`;
  viewer.classList.add('empty');
  $(which === 'seq1' ? 'summary1' : 'summary2').textContent = 'No analysis yet.';
  const bank = $(which === 'seq1' ? 'bank1' : 'bank2');
  bank.className = 'bank empty-bank';
  bank.textContent = `Digest ${which === 'seq1' ? 'sequence 1' : 'sequence 2'} to generate draggable fragments.`;
  state.analyses[which] = null;
}

function copyAssembledSequence() {
  if (!state.assembly.length) return;
  const mergedSeq = state.assembly.map(item => item.variant.seq).join('');
  navigator.clipboard.writeText(mergedSeq);
  $('assemblyStatus').textContent = 'Assembled sequence copied to clipboard.';
}

function init() {
  populateEnzymeSelect($('enzymes1'));
  populateEnzymeSelect($('enzymes2'));
  populateEnzymeSelect($('builderEnzyme'));
  ['BsaI', 'BsmBI', 'SapI'].forEach(name => {
    for (const selectId of ['enzymes1', 'enzymes2']) {
      const opt = Array.from($(selectId).options).find(o => o.value === name);
      if (opt) opt.selected = true;
    }
  });
  $('builderEnzyme').value = 'BsaI';
  readFileIntoInput('seq1File', 'seq1Input');
  readFileIntoInput('seq2File', 'seq2Input');
  $('analyze1').addEventListener('click', () => analyzeSequence('seq1'));
  $('analyze2').addEventListener('click', () => analyzeSequence('seq2'));
  $('digest1').addEventListener('click', () => digestSequence('seq1'));
  $('digest2').addEventListener('click', () => digestSequence('seq2'));
  $('clear1').addEventListener('click', () => clearSequence('seq1'));
  $('clear2').addEventListener('click', () => clearSequence('seq2'));
  $('buildSite').addEventListener('click', applyBuilder);
  $('copyBlock').addEventListener('click', copyBlockOnly);
  $('clearWorkbench').addEventListener('click', () => { state.assembly = []; renderAssemblyTrack(); });
  $('copyAssembled').addEventListener('click', copyAssembledSequence);
  renderSyntheticBank();
  renderAssemblyTrack();
}

window.addEventListener('DOMContentLoaded', init);
