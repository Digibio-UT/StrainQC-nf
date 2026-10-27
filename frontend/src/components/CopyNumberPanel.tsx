import React from 'react';
import {
    BarChart,
    Bar,
    XAxis,
    YAxis,
    CartesianGrid,
    Tooltip,
    ResponsiveContainer,
    Cell,
    ReferenceLine
} from 'recharts';

// Thresholds mirror bin/summarize_cnv.py (chromosome mean CN relative to the genome median).
const GAIN_RATIO = 1.20;
const LOSS_RATIO = 0.75;

export interface ChromosomeCopyNumber {
    seq_id: string;
    mean_cn: number;
    cn_vs_genome_median: number;
    length: number;
    n_segments: number;
    altered_fraction: number;
    call: 'gain' | 'loss' | 'focal' | 'neutral';
    aneuploidy: boolean;
}

interface CopyNumberPanelProps {
    data: ChromosomeCopyNumber[];
    aneuploidy: { seq_id: string; call: string; mean_cn: number }[];
    tool?: string;
    displayName: (id: string) => string;
    selectedChrom: string;
    onSelectChrom: (seqId: string) => void;
}

const CALL_COLORS: Record<string, string> = {
    gain: '#d32f2f',
    loss: '#1565c0',
    focal: '#ffcc80',   // muted: sub-chromosomal segments are context, not the headline
    neutral: '#90a4ae',
};

const CALL_LABELS: Record<string, string> = {
    gain: 'Whole-chromosome gain',
    loss: 'Whole-chromosome loss',
    focal: 'Partial change (segment)',
    neutral: 'No change',
};

const CopyNumberPanel: React.FC<CopyNumberPanelProps> = ({
    data, aneuploidy, tool, displayName: accessionName, selectedChrom, onSelectChrom
}) => {
    // Ensembl-style references name chromosomes plain "I", "II", ...; show them as "Chr I".
    const displayName = (id: string) => {
        const n = accessionName(id);
        return /^[IVX]+$/.test(n) ? `Chr ${n}` : n;
    };

    const chartData = data.map(c => ({
        ...c,
        label: displayName(c.seq_id).replace(/^Chr /, ''),
        color:CALL_COLORS[c.call] || CALL_COLORS.neutral,
    }));

    const ratioOf = (seqId: string) => data.find(c => c.seq_id === seqId)?.cn_vs_genome_median;

    const yMax =Math.max(1.5, ...data.map(c => c.cn_vs_genome_median + 0.1));

    const renderTooltip = ({ active, payload }: { active?: boolean; payload?: readonly { payload?: typeof chartData[number] }[] }) => {
        const c = active && payload && payload.length ? payload[0].payload : undefined;
        if (!c) return null;
        return (
            <div style={{ backgroundColor: '#fff', padding: '10px', border: '1px solid #ccc', borderRadius: '4px', fontSize: '0.85em' }}>
                <p style={{ margin: 0, fontWeight: 'bold' }}>
                    {displayName(c.seq_id)}{accessionName(c.seq_id) !== c.seq_id && ` (${c.seq_id})`}
                </p>
                <p style={{ margin: '4px 0 0 0', color: c.color, fontWeight: 600 }}>{CALL_LABELS[c.call] || c.call}</p>
                <p style={{ margin: '4px 0 0 0' }}>Relative copy number: {c.cn_vs_genome_median.toFixed(2)}×</p>
                <p style={{ margin: 0, color: '#666' }}>Estimated copies: {c.mean_cn.toFixed(2)}</p>
                <p style={{ margin: 0, color: '#666' }}>Segments: {c.n_segments} · altered: {Math.round(c.altered_fraction * 100)}% of length</p>
                <p style={{ margin: 0, color: '#666' }}>Length: {c.length.toLocaleString()} bp</p>
            </div>
        );
    };

    return (
        <div style={{ marginBottom: '40px', padding: '20px', backgroundColor: '#fff', borderRadius: '8px', boxShadow: '0 2px 4px rgba(0,0,0,0.1)' }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: '12px', flexWrap: 'wrap' }}>
                <h3 style={{ color: '#2c3e50', margin: 0 }}>Copy Number &amp; Aneuploidy</h3>
                <span style={{
                    fontSize: '0.75em',
                    padding: '3px 8px',
                    borderRadius: '12px',
                    backgroundColor: '#f3e5f5',
                    color: '#6a1b9a',
                    border: '1px solid #ce93d8',
                    fontWeight: 600,
                }}>
                    {tool === 'cnvkit' || !tool ? 'CNVkit' : tool} · parent as reference
                </span>
            </div>

            {aneuploidy.length > 0 ? (
                <div style={{ marginTop: '12px', padding: '10px 14px', backgroundColor: '#ffebee', border: '1px solid #ef9a9a', borderRadius: '6px' }}>
                    <strong style={{ color: '#c62828' }}>Aneuploidy detected: </strong>
                    {aneuploidy.map((a, i) => (
                        <span key={a.seq_id}>
                            {i > 0 && ', '}
                            <span
                                onClick={() => onSelectChrom(a.seq_id)}
                                style={{ cursor: 'pointer', textDecoration: 'underline' }}
                                title="Show this chromosome in the coverage plot"
                            >
                                {displayName(a.seq_id)} {a.call}
                            </span>
                            {ratioOf(a.seq_id) != null && ` (${ratioOf(a.seq_id)!.toFixed(2)}× the rest of the genome)`}
                        </span>
                    ))}
                </div>
            ) : (
                <div style={{ marginTop: '12px', padding: '10px 14px', backgroundColor: '#e8f5e9', border: '1px solid #a5d6a7', borderRadius: '6px' }}>
                    <strong style={{ color: '#2e7d32' }}>No whole-chromosome gains or losses</strong>
                    <span style={{ color: '#555' }}> compared with the parent strain.</span>
                </div>
            )}

            <p style={{ fontSize: '0.85em', color: '#888', margin: '12px 0 10px 0' }}>
                Copy number of each chromosome relative to the genome as a whole (1.0 = same as the parent).
                Bars above {GAIN_RATIO.toFixed(2)}× or below {LOSS_RATIO.toFixed(2)}× over most of a chromosome are called as aneuploidy.
                Click a bar to open that chromosome in the coverage plot.
            </p>

            <div className="cnv-chart" style={{ width: '100%', height: 260 }}>
                <style>{`.cnv-chart *:focus { outline: none; }`}</style>
                <ResponsiveContainer width="100%" height="100%">
                    <BarChart data={chartData} margin={{ top: 10, right: 30, left: 10, bottom: 5 }}>
                        <CartesianGrid strokeDasharray="3 3" vertical={false} />
                        <XAxis dataKey="label" tick={{ fontSize: 12 }} interval={0} />
                        <YAxis
                            domain={[0, Math.ceil(yMax * 4) / 4]}
                            tickFormatter={(v: number) => `${v}×`}
                            tick={{ fontSize: 12 }}
                            label={{ value: 'Relative copy number', angle: -90, position: 'insideLeft', style: { fontSize: 12, fill: '#666', textAnchor: 'middle' } }}
                        />
                        <Tooltip content={renderTooltip} cursor={{ fill: 'rgba(0,0,0,0.05)' }} />
                        <ReferenceLine y={1} stroke="#555" />
                        <ReferenceLine y={GAIN_RATIO} stroke={CALL_COLORS.gain} strokeDasharray="4 4"
                            label={{ value: 'gain', position: 'right', fill: CALL_COLORS.gain, fontSize: 11 }} />
                        <ReferenceLine y={LOSS_RATIO} stroke={CALL_COLORS.loss} strokeDasharray="4 4"
                            label={{ value: 'loss', position: 'right', fill: CALL_COLORS.loss, fontSize: 11 }} />
                        <Bar
                            dataKey="cn_vs_genome_median"
                            isAnimationActive={false}
                            cursor="pointer"
                            onClick={(d) => {
                                const id = d?.payload?.seq_id;
                                if (id) onSelectChrom(id);
                            }}
                        >
                            {chartData.map(c => (
                                <Cell
                                    key={c.seq_id}
                                    fill={c.color}
                                    stroke={c.seq_id === selectedChrom ? '#212121' : 'none'}
                                    strokeWidth={c.seq_id === selectedChrom ? 2 : 0}
                                />
                            ))}
                        </Bar>
                    </BarChart>
                </ResponsiveContainer>
            </div>

            <div style={{ display: 'flex', gap: '16px', flexWrap: 'wrap', fontSize: '0.8em', color: '#555', marginTop: '6px' }}>
                {(['gain', 'loss', 'focal', 'neutral'] as const).map(k => (
                    <span key={k} style={{ display: 'flex', alignItems: 'center', gap: '5px' }}>
                        <span style={{ width: '12px', height: '12px', backgroundColor: CALL_COLORS[k], borderRadius: '2px', display: 'inline-block' }} />
                        {CALL_LABELS[k]}
                    </span>
                ))}
            </div>
        </div>
    );
};

export default CopyNumberPanel;
