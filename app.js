const IUPAC = {
  A: 'A', C: 'C', G: 'G', T: 'T', U: 'T',
  R: 'AG', Y: 'CT', S: 'GC', W: 'AT', K: 'GT', M: 'AC',
  B: 'CGT', D: 'AGT', H: 'ACT', V: 'ACG', N: 'ACGT'
};

const state = {
  analyses: {
    seq1: null,
    seq2: null
  }
};

function $(id) { return document.getElementById(id); }

function revComp(seq) {
  const map = { A: 'T', T: 'A', G: 'C', C: 'G', R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K', B: 'V', D: 'H', H: 'D', V: 'B', N: 'N' };
  return seq.toUpperCase().split('').reverse().map(ch => map[ch] || 'N').join('');
}

function iupacToRegex(site) {
  return site.toUpperCase().split('').map(ch => `[${IUPAC[ch] || 'ACGT'}]`).join('');
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

  return text
    .split(/\r?\n/)
    .filter(line => !line.startsWith('>'))
    .join('')
    .toUpperCase()
    .replace(/U/g, 'T')
    .replace(/[^ACGTRYSWKMBDHVN]/g, '');
}

function formatCut(enzyme) {
  const c = enzyme.cut;
  if (c.kind === 'simple') return `${enzyme.recognition}(${c.top}/${c.bottom})`;
  return `(${c.leftTop}/${c.leftBottom})${enzyme.recognition}(${c.rightTop}/${c.rightBottom})`;
}

function populateEnzymeSelect(select) {
  select.innerHTML = '';
  window.TYPE_IIS_ENZYMES
    .slice()
    .sort((a, b) => a.name.localeCompare(b.name))
    .forEach(enzyme => {
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
    for (let i = 0; i <= upper.length - L; i += 1) {
      const fragment = upper.slice(i, i + L);
      let forward = true;
      for (let j = 0; j < L; j += 1) {
        if (!baseMatches(fragment[j], site[j])) { forward = false; break; }
      }
      if (forward) hits.push(buildHit(upper, enzyme, i, 'forward'));

      const rc = revComp(site);
      let reverse = true;
      for (let j = 0; j < L; j += 1) {
        if (!baseMatches(fragment[j], rc[j])) { reverse = false; break; }
      }
      if (reverse) hits.push(buildHit(upper, enzyme, i, 'reverse'));
    }
  }

  return hits.sort((a, b) => a.start - b.start || a.enzyme.name.localeCompare(b.enzyme.name) || a.orientation.localeCompare(b.orientation));
}

function buildHit(seq, enzyme, start, orientation) {
  const L = enzyme.recognition.length;
  const end = start + L;
  const hit = {
    enzyme,
    orientation,
    start,
    end,
    recognitionShown: orientation === 'forward' ? enzyme.recognition : revComp(enzyme.recognition),
    cuts: []
  };

  if (enzyme.cut.kind === 'simple') {
    let topCut, bottomCut;
    if (orientation === 'forward') {
      topCut = end + enzyme.cut.top;
      bottomCut = end + enzyme.cut.bottom;
    } else {
      topCut = start - enzyme.cut.bottom;
      bottomCut = start - enzyme.cut.top;
    }
    hit.cuts.push(makeCutDescriptor(seq, topCut, bottomCut, enzyme, orientation, start, end));
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
    hit.cuts.push(makeCutDescriptor(seq, leftTop, leftBottom, enzyme, orientation, start, end, 'left'));
    hit.cuts.push(makeCutDescriptor(seq, rightTop, rightBottom, enzyme, orientation, start, end, 'right'));
  }

  hit.display = renderHitDisplay(seq, hit);
  return hit;
}

function makeCutDescriptor(seq, topCut, bottomCut, enzyme, orientation, start, end, flank = 'single') {
  const leftEnd = describeFragmentEnd(seq, topCut, bottomCut, 'left');
  const rightEnd = describeFragmentEnd(seq, topCut, bottomCut, 'right');
  return {
    flank,
    topCut,
    bottomCut,
    inBounds: topCut >= 0 && bottomCut >= 0 && topCut <= seq.length && bottomCut <= seq.length,
    overhangLength: Math.abs(topCut - bottomCut),
    leftEnd,
    rightEnd,
    id: `${enzyme.name}|${orientation}|${start}|${flank}|${topCut}|${bottomCut}`
  };
}

function describeFragmentEnd(seq, topCut, bottomCut, side) {
  const diff = topCut - bottomCut;
  if (diff === 0) {
    return { side, type: 'blunt', length: 0, sequence: '', compatibleKey: 'blunt' };
  }

  if (side === 'left') {
    if (diff > 0) {
      const over = seq.slice(bottomCut, topCut);
      return { side, type: '5prime', length: diff, sequence: over, compatibleKey: `5:${over}` };
    }
    const over = revComp(seq.slice(topCut, bottomCut));
    return { side, type: '3prime', length: -diff, sequence: over, compatibleKey: `3:${over}` };
  }

  if (diff < 0) {
    const over = seq.slice(topCut, bottomCut);
    return { side, type: '5prime', length: -diff, sequence: over, compatibleKey: `5:${over}` };
  }
  const over = revComp(seq.slice(bottomCut, topCut));
  return { side, type: '3prime', length: diff, sequence: over, compatibleKey: `3:${over}` };
}

function areEndsCompatible(a, b) {
  if (!a || !b) return { ok: false, reason: 'Missing end selection.' };
  if (a.type === 'blunt' && b.type === 'blunt') return { ok: true, reason: 'Both ends are blunt.' };
  if (a.type !== b.type) return { ok: false, reason: 'One end is 5′ overhang and the other is 3′ overhang.' };
  if (a.length !== b.length) return { ok: false, reason: 'Overhang lengths differ.' };
  const comp = revComp(b.sequence);
  if (a.sequence === comp) return { ok: true, reason: `Complementary ${a.type === '5prime' ? '5′' : '3′'} overhangs.` };
  return { ok: false, reason: 'Overhang sequences are not reverse complements.' };
}

function renderHitDisplay(seq, hit) {
  const padLeft = Math.max(0, hit.start - 8);
  const padRight = Math.min(seq.length, hit.end + 24);
  const top = seq.slice(padLeft, padRight);
  const bottom = revComp(top);
  const lines = [];
  lines.push(`Context ${padLeft + 1}-${padRight}`);
  for (const cut of hit.cuts) {
    const topMark = buildMarkerLine(padLeft, padRight, cut.topCut);
    const botMark = buildMarkerLine(padLeft, padRight, cut.bottomCut);
    lines.push(`5' ${top}`);
    lines.push(`   ${topMark}  top cut`);
    lines.push(`3' ${bottom}`);
    lines.push(`   ${botMark}  bottom cut`);
    lines.push(`Cut ${cut.flank}: top ${cut.topCut}, bottom ${cut.bottomCut} | left end ${describeEnd(cut.leftEnd)} | right end ${describeEnd(cut.rightEnd)}`);
    lines.push('');
  }
  return lines.join('\n');
}

function buildMarkerLine(windowStart, windowEnd, cutPos) {
  const width = windowEnd - windowStart;
  const chars = new Array(width).fill(' ');
  const idx = cutPos - windowStart;
  if (idx >= 0 && idx < width) chars[idx] = '│';
  return chars.join('');
}

function describeEnd(end) {
  if (end.type === 'blunt') return 'blunt';
  return `${end.type === '5prime' ? '5′' : '3′'} ${end.length} nt (${end.sequence || '∅'})`;
}

function analyzeSequence(which) {
  const seqRaw = $(`${which}Input`).value;
  const name = $(`${which}Name`).value.trim() || which;
  const enzymes = getSelectedEnzymes(which === 'seq1' ? 'enzymes1' : 'enzymes2');
  const seq = cleanSequence(seqRaw);

  if (!seq) {
    alert(`No valid sequence found in ${name}.`);
    return;
  }
  if (!enzymes.length) {
    alert('Select at least one Type IIS enzyme first.');
    return;
  }

  const hits = findSites(seq, enzymes);
  const simpleEnds = [];
  hits.forEach(hit => hit.cuts.forEach(c => {
    if (hit.enzyme.cut.kind === 'simple' && c.inBounds) {
      simpleEnds.push({
        id: `${which}:${c.id}:left`,
        label: `${hit.enzyme.name} ${hit.orientation} site ${hit.start + 1}-${hit.end} | left fragment end | ${describeEnd(c.leftEnd)}`,
        end: c.leftEnd,
        meta: { enzyme: hit.enzyme.name, orientation: hit.orientation, siteStart: hit.start + 1, siteEnd: hit.end, flank: c.flank }
      });
      simpleEnds.push({
        id: `${which}:${c.id}:right`,
        label: `${hit.enzyme.name} ${hit.orientation} site ${hit.start + 1}-${hit.end} | right fragment end | ${describeEnd(c.rightEnd)}`,
        end: c.rightEnd,
        meta: { enzyme: hit.enzyme.name, orientation: hit.orientation, siteStart: hit.start + 1, siteEnd: hit.end, flank: c.flank }
      });
    }
  }));

  state.analyses[which] = { which, name, seq, hits, simpleEnds };
  renderAnalysis(which);
  refreshEndSelectors();
  renderGlobalSummary();
}

function renderAnalysis(which) {
  const targetSummary = $(which === 'seq1' ? 'summary1' : 'summary2');
  const targetHits = $(which === 'seq1' ? 'hits1' : 'hits2');
  const analysis = state.analyses[which];
  if (!analysis) return;

  targetSummary.innerHTML = `${analysis.name}: ${analysis.seq.length.toLocaleString()} nt, ${analysis.hits.length} detected site(s), ${analysis.simpleEnds.length} simple ligatable end(s).`;
  targetHits.innerHTML = '';

  if (!analysis.hits.length) {
    targetHits.innerHTML = '<div class="hit">No selected Type IIS sites were found in this sequence.</div>';
    return;
  }

  for (const hit of analysis.hits) {
    const div = document.createElement('div');
    div.className = 'hit';
    const cutNotes = hit.cuts.map(c => `${c.flank}: top ${c.topCut}, bottom ${c.bottomCut}${c.inBounds ? '' : ' (out of bounds for this sequence window)'}`).join(' · ');
    div.innerHTML = `
      <div class="row"><div><strong>${hit.enzyme.name}</strong> <span class="badge">${hit.orientation}</span> <span class="badge">${formatCut(hit.enzyme)}</span></div></div>
      <div class="kv">
        <div>Recognition</div><div>${hit.recognitionShown}</div>
        <div>Sequence span</div><div>${hit.start + 1}-${hit.end}</div>
        <div>Cut positions</div><div>${cutNotes}</div>
        <div>Support</div><div>${hit.enzyme.cut.kind === 'simple' ? 'Full digest + ligation support' : 'Detection and cut-display support; dual-flank compatibility not enumerated'}</div>
      </div>
      <pre>${escapeHtml(hit.display)}</pre>
    `;
    targetHits.appendChild(div);
  }
}

function escapeHtml(text) {
  return text.replace(/[&<>]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' }[c]));
}

function refreshEndSelectors() {
  for (const id of ['end1', 'end2']) $(id).innerHTML = '';
  const map = { end1: state.analyses.seq1, end2: state.analyses.seq2 };

  for (const [id, analysis] of Object.entries(map)) {
    const select = $(id);
    const placeholder = document.createElement('option');
    placeholder.value = '';
    placeholder.textContent = analysis ? `Choose an end from ${analysis.name}` : 'Analyze the sequence first';
    select.appendChild(placeholder);
    if (!analysis) continue;
    analysis.simpleEnds.forEach(item => {
      const opt = document.createElement('option');
      opt.value = item.id;
      opt.textContent = item.label;
      select.appendChild(opt);
    });
  }
}

function findEndById(id) {
  for (const key of ['seq1', 'seq2']) {
    const analysis = state.analyses[key];
    const found = analysis?.simpleEnds.find(item => item.id === id);
    if (found) return found;
  }
  return null;
}

function compareSelectedEnds() {
  const a = findEndById($('end1').value);
  const b = findEndById($('end2').value);
  const result = areEndsCompatible(a?.end, b?.end);
  $('compatibilityResult').innerHTML = `
    <div class="status ${result.ok ? 'ok' : 'no'}">${result.ok ? 'Compatible' : 'Not compatible'}</div>
    <div style="margin-top:8px;">${result.reason}</div>
    <div class="small" style="margin-top:8px;">${a ? a.label : 'No end selected from sequence 1'}<br>${b ? b.label : 'No end selected from sequence 2'}</div>
  `;
}

function scanAllPairs() {
  const a = state.analyses.seq1?.simpleEnds || [];
  const b = state.analyses.seq2?.simpleEnds || [];
  if (!a.length || !b.length) {
    $('compatibilityResult').innerHTML = '<div class="status warn">Analyze both sequences first and make sure they produce simple ends within bounds.</div>';
    return;
  }
  const matches = [];
  for (const left of a) {
    for (const right of b) {
      const res = areEndsCompatible(left.end, right.end);
      if (res.ok) matches.push({ left, right, reason: res.reason });
    }
  }
  if (!matches.length) {
    $('compatibilityResult').innerHTML = '<div class="status no">No compatible ligation pairs were found between the enumerated simple ends.</div>';
    return;
  }
  $('compatibilityResult').innerHTML = `
    <div class="status ok">${matches.length} compatible pair(s) found</div>
    <div class="small" style="margin-top:8px;">${matches.slice(0, 20).map(m => `• ${m.left.label}<br>&nbsp;&nbsp;&nbsp;↔ ${m.right.label}`).join('<br><br>')}</div>
  `;
}

function renderGlobalSummary() {
  const a = state.analyses.seq1;
  const b = state.analyses.seq2;
  const lines = [];
  lines.push(`${window.TYPE_IIS_ENZYMES.length} curated Type IIS enzyme entries are loaded.`);
  lines.push('Both strand orientations are scanned for each selected recognition site.');
  lines.push('Simple cutters are supported for end-compatibility calculations; dual-flank cutters are displayed and mapped but not enumerated into ligation pairs.');
  if (a) lines.push(`${a.name}: ${a.hits.length} site(s), ${a.simpleEnds.length} simple end(s).`);
  if (b) lines.push(`${b.name}: ${b.hits.length} site(s), ${b.simpleEnds.length} simple end(s).`);
  $('globalSummary').innerHTML = lines.join('<br>');
}

function buildRestrictionBlock() {
  const enzyme = window.TYPE_IIS_ENZYMES.find(e => e.name === $('builderEnzyme').value);
  const left = cleanSequence($('builderLeft').value);
  const payload = cleanSequence($('builderPayload').value);
  const right = cleanSequence($('builderRight').value);
  const block = `${left}${enzyme.recognition}${payload}${right}`;
  $('builderPreview').innerHTML = `Block: <code>${block}</code><br><span class="small">Pattern: ${formatCut(enzyme)}</span>`;
  return block;
}

function applyBlock() {
  const block = buildRestrictionBlock();
  const mode = $('builderMode').value;
  const current = cleanSequence($('seq1Input').value);
  if (mode === 'replace') $('seq1Input').value = block;
  if (mode === 'prefix') $('seq1Input').value = block + current;
  if (mode === 'suffix') $('seq1Input').value = current + block;
}

async function copyBlockOnly() {
  const block = buildRestrictionBlock();
  try {
    await navigator.clipboard.writeText(block);
    $('builderPreview').innerHTML += '<br>Copied to clipboard.';
  } catch {
    $('builderPreview').innerHTML += '<br>Clipboard copy failed in this browser context.';
  }
}

function loadFile(inputId, targetTextareaId) {
  const file = $(inputId).files?.[0];
  if (!file) return;
  const reader = new FileReader();
  reader.onload = () => { $(targetTextareaId).value = String(reader.result || ''); };
  reader.readAsText(file);
}

function init() {
  ['enzymes1', 'enzymes2', 'builderEnzyme'].forEach(id => populateEnzymeSelect($(id)));
  ['BsaI', 'BsmBI', 'SapI'].forEach(name => {
    const opt1 = Array.from($('enzymes1').options).find(o => o.value === name);
    const opt2 = Array.from($('enzymes2').options).find(o => o.value === name);
    if (opt1) opt1.selected = true;
    if (opt2) opt2.selected = true;
  });
  $('builderEnzyme').value = 'BsaI';
  buildRestrictionBlock();
  renderGlobalSummary();

  $('analyze1').addEventListener('click', () => analyzeSequence('seq1'));
  $('analyze2').addEventListener('click', () => analyzeSequence('seq2'));
  $('compareEnds').addEventListener('click', compareSelectedEnds);
  $('scanPairs').addEventListener('click', scanAllPairs);
  $('buildSite').addEventListener('click', applyBlock);
  $('copyBlock').addEventListener('click', copyBlockOnly);
  $('builderEnzyme').addEventListener('change', buildRestrictionBlock);
  ['builderLeft', 'builderPayload', 'builderRight'].forEach(id => $(id).addEventListener('input', buildRestrictionBlock));
  $('seq1File').addEventListener('change', () => loadFile('seq1File', 'seq1Input'));
  $('seq2File').addEventListener('change', () => loadFile('seq2File', 'seq2Input'));
  $('clear1').addEventListener('click', () => { $('seq1Input').value = ''; $('summary1').textContent = 'No analysis yet.'; $('hits1').innerHTML = ''; state.analyses.seq1 = null; refreshEndSelectors(); renderGlobalSummary(); });
  $('clear2').addEventListener('click', () => { $('seq2Input').value = ''; $('summary2').textContent = 'No analysis yet.'; $('hits2').innerHTML = ''; state.analyses.seq2 = null; refreshEndSelectors(); renderGlobalSummary(); });
}

init();
