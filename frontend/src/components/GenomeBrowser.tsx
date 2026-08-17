import React, { useRef, useState, useCallback, useMemo, useEffect } from 'react';

// ─── Types ───────────────────────────────────────────────────────────────────

export interface GenomeFeature {
    start: number;
    end: number;
    strand: number;
    name: string;
    locus_tag: string;
    product: string;
    seq_id?: string;
    original_id?: string;
    parts?: { start: number; end: number }[];
}

export interface VariantMarker {
    pos: number;
    end: number;
    size: number;
    type: string;
    mutation: string;
    gene: string;
    annotation: string;
    freq: number;
    fill: string;
    evidence_links?: { text: string; url: string }[];
}

interface GenomeBrowserProps {
    features: GenomeFeature[];
    variants: VariantMarker[];
    viewStart: number;
    viewEnd: number;
    maxPos: number;
    onViewChange: (start: number, end: number) => void;
    onFeatureClick?: (feature: GenomeFeature) => void;
    onVariantClick?: (variant: VariantMarker) => void;
    selectedVariantPos?: number | null;
    jobId?: string;
    onEvidenceOpen?: (url: string) => void;
}

// ─── Constants ───────────────────────────────────────────────────────────────

const MARGIN = { top: 10, right: 20, bottom: 30, left: 50 };
const GENE_HEIGHT = 16;
const GENE_GAP = 2;
const GENE_ROW_STEP = GENE_HEIGHT + GENE_GAP;
const MIN_STRAND_SPACE = 40; // minimum pixels for each strand section
const VARIANT_BAND = 30;     // pixels reserved for variant markers
const ARROW_SIZE = 8;

// ─── Helpers ─────────────────────────────────────────────────────────────────

function formatBp(value: number): string {
    if (value === 0) return '0';
    if (Math.abs(value) >= 1_000_000) return `${(value / 1_000_000).toFixed(2)}M`;
    if (Math.abs(value) >= 1_000) return `${(value / 1_000).toFixed(1)}k`;
    return value.toString();
}

function generateTicks(start: number, end: number, maxTicks: number = 8): number[] {
    const range = end - start;
    if (range <= 0) return [start];
    
    // Find a nice tick interval
    const rawInterval = range / maxTicks;
    const magnitude = Math.pow(10, Math.floor(Math.log10(rawInterval)));
    const residual = rawInterval / magnitude;
    
    let interval: number;
    if (residual <= 1.5) interval = magnitude;
    else if (residual <= 3.5) interval = 2 * magnitude;
    else if (residual <= 7.5) interval = 5 * magnitude;
    else interval = 10 * magnitude;
    
    const ticks: number[] = [];
    const first = Math.ceil(start / interval) * interval;
    for (let t = first; t <= end; t += interval) {
        ticks.push(t);
    }
    return ticks;
}

// ─── Component ───────────────────────────────────────────────────────────────

const GenomeBrowser: React.FC<GenomeBrowserProps> = ({
    features,
    variants,
    viewStart,
    viewEnd,
    maxPos,
    onViewChange,
    onFeatureClick,
    onVariantClick,
    selectedVariantPos,
    onEvidenceOpen,
}) => {
    const svgRef = useRef<SVGSVGElement>(null);
    const containerRef = useRef<HTMLDivElement>(null);
    const [svgWidth, setSvgWidth] = useState(800);
    const [hoveredFeature, setHoveredFeature] = useState<GenomeFeature | null>(null);
    const [hoveredVariant, setHoveredVariant] = useState<VariantMarker | null>(null);
    const [tooltipPos, setTooltipPos] = useState({ x: 0, y: 0 });
    const [isDragging, setIsDragging] = useState(false);
    const dragStartRef = useRef<{ x: number; viewStart: number; viewEnd: number } | null>(null);

    // Responsive width
    useEffect(() => {
        const measure = () => {
            if (containerRef.current) {
                setSvgWidth(containerRef.current.clientWidth);
            }
        };
        measure();
        const observer = new ResizeObserver(measure);
        if (containerRef.current) observer.observe(containerRef.current);
        return () => observer.disconnect();
    }, []);

    const plotWidth = svgWidth - MARGIN.left - MARGIN.right;
    const viewSpan = viewEnd - viewStart;

    // Scale: bp -> pixel x
    const bpToX = useCallback((bp: number) => {
        return MARGIN.left + ((bp - viewStart) / (viewSpan || 1)) * plotWidth;
    }, [viewStart, viewSpan, plotWidth]);

    // Scale: pixel x -> bp
    const xToBp = useCallback((px: number) => {
        return viewStart + ((px - MARGIN.left) / plotWidth) * viewSpan;
    }, [viewStart, viewSpan, plotWidth]);

    // ── Filtered data ────────────────────────────────────────────────────────
    const buffer = viewSpan * 0.1;
    const visibleFeatures = useMemo(() => {
        return features.filter(f =>
            f.end > viewStart - buffer && f.start < viewEnd + buffer
        );
    }, [features, viewStart, viewEnd, buffer]);

    const visibleVariants = useMemo(() => {
        return variants.filter(v =>
            v.pos >= viewStart - buffer && v.pos <= viewEnd + buffer
        );
    }, [variants, viewStart, viewEnd, buffer]);

    const visibleDeletions = useMemo(() => {
        return variants.filter(v =>
            v.type.includes('DEL') &&
            (v.end || v.pos) > viewStart && v.pos < viewEnd
        );
    }, [variants, viewStart, viewEnd]);

    // ── Layout genes to avoid overlap ────────────────────────────────────────
    const assignGeneRows = useCallback((genes: GenomeFeature[]) => {
        // Simple row-based stacking — returns row assignments without Y positions
        const rows: { end: number }[] = [];
        const laid: { gene: GenomeFeature; row: number; x1: number; x2: number }[] = [];

        const sorted = [...genes].sort((a, b) => a.start - b.start);
        for (const gene of sorted) {
            const x1 = bpToX(Math.max(gene.start, viewStart));
            const x2 = bpToX(Math.min(gene.end, viewEnd));
            if (x2 - x1 < 0.5) continue; // too small to render

            // Find first row where this gene fits
            let row = 0;
            for (let r = 0; r < rows.length; r++) {
                if (rows[r].end + 4 < x1) {  // 4px gap between genes
                    row = r;
                    break;
                }
                row = r + 1;
            }
            if (row >= rows.length) rows.push({ end: x2 });
            else rows[row].end = x2;

            laid.push({ gene, row, x1, x2 });
        }
        return { items: laid, rowCount: rows.length };
    }, [bpToX, viewStart, viewEnd]);

    const forwardLayout = useMemo(() =>
        assignGeneRows(visibleFeatures.filter(f => f.strand === 1)),
    [visibleFeatures, assignGeneRows]);

    const reverseLayout = useMemo(() =>
        assignGeneRows(visibleFeatures.filter(f => f.strand === -1)),
    [visibleFeatures, assignGeneRows]);

    // ── Dynamic track sizing ─────────────────────────────────────────────────
    const forwardSpace = Math.max(MIN_STRAND_SPACE, forwardLayout.rowCount * GENE_ROW_STEP + 20);
    const reverseSpace = Math.max(MIN_STRAND_SPACE, reverseLayout.rowCount * GENE_ROW_STEP + 20);
    const trackHeight = forwardSpace + VARIANT_BAND + reverseSpace;

    // Vertical anchors
    const FORWARD_Y = forwardSpace - 10;       // bottom of forward section (genes stack upward)
    const VARIANT_Y = forwardSpace + VARIANT_BAND / 2;
    const REVERSE_Y = forwardSpace + VARIANT_BAND + 10; // top of reverse section (genes stack downward)

    // Assign final Y positions
    const forwardGenes = useMemo(() =>
        forwardLayout.items.map(item => ({
            ...item,
            y: FORWARD_Y - item.row * GENE_ROW_STEP,
        })),
    [forwardLayout, FORWARD_Y]);

    const reverseGenes = useMemo(() =>
        reverseLayout.items.map(item => ({
            ...item,
            y: REVERSE_Y + item.row * GENE_ROW_STEP,
        })),
    [reverseLayout, REVERSE_Y]);

    // ── Interaction handlers ─────────────────────────────────────────────────

    const handleWheel = useCallback((e: React.WheelEvent) => {
        e.preventDefault();
        const zoomFactor = e.deltaY > 0 ? 1.2 : 1 / 1.2;
        
        // Zoom centered on mouse position
        const rect = svgRef.current?.getBoundingClientRect();
        if (!rect) return;
        const mouseX = e.clientX - rect.left;
        const mouseBp = xToBp(mouseX);

        const newSpan = Math.max(200, Math.min(maxPos, viewSpan * zoomFactor));
        const ratio = (mouseBp - viewStart) / (viewSpan || 1);
        let newStart = mouseBp - ratio * newSpan;
        let newEnd = mouseBp + (1 - ratio) * newSpan;

        // Clamp
        if (newStart < 0) { newEnd -= newStart; newStart = 0; }
        if (newEnd > maxPos) { newStart -= (newEnd - maxPos); newEnd = maxPos; }
        newStart = Math.max(0, newStart);

        onViewChange(Math.round(newStart), Math.round(newEnd));
    }, [viewStart, viewSpan, maxPos, xToBp, onViewChange]);

    const handleMouseDown = useCallback((e: React.MouseEvent) => {
        if (e.button !== 0) return; // left click only
        setIsDragging(true);
        dragStartRef.current = {
            x: e.clientX,
            viewStart,
            viewEnd,
        };
    }, [viewStart, viewEnd]);

    const handleMouseMove = useCallback((e: React.MouseEvent) => {
        if (!isDragging || !dragStartRef.current) return;
        const dx = e.clientX - dragStartRef.current.x;
        const bpDelta = -(dx / plotWidth) * viewSpan;

        let newStart = dragStartRef.current.viewStart + bpDelta;
        let newEnd = dragStartRef.current.viewEnd + bpDelta;

        // Clamp
        if (newStart < 0) { newEnd -= newStart; newStart = 0; }
        if (newEnd > maxPos) { newStart -= (newEnd - maxPos); newEnd = maxPos; }
        newStart = Math.max(0, newStart);
        newEnd = Math.min(maxPos, newEnd);

        onViewChange(Math.round(newStart), Math.round(newEnd));
    }, [isDragging, plotWidth, viewSpan, maxPos, onViewChange]);

    const handleMouseUp = useCallback(() => {
        setIsDragging(false);
        dragStartRef.current = null;
    }, []);

    // Release drag if mouse leaves
    useEffect(() => {
        const handleGlobalUp = () => {
            setIsDragging(false);
            dragStartRef.current = null;
        };
        window.addEventListener('mouseup', handleGlobalUp);
        return () => window.removeEventListener('mouseup', handleGlobalUp);
    }, []);

    // ── Tooltip handlers ─────────────────────────────────────────────────────
    const handleFeatureHover = useCallback((gene: GenomeFeature, e: React.MouseEvent) => {
        setHoveredFeature(gene);
        setHoveredVariant(null);
        const rect = containerRef.current?.getBoundingClientRect();
        if (rect) {
            setTooltipPos({ x: e.clientX - rect.left, y: e.clientY - rect.top });
        }
    }, []);

    const handleVariantHover = useCallback((variant: VariantMarker, e: React.MouseEvent) => {
        setHoveredVariant(variant);
        setHoveredFeature(null);
        const rect = containerRef.current?.getBoundingClientRect();
        if (rect) {
            setTooltipPos({ x: e.clientX - rect.left, y: e.clientY - rect.top });
        }
    }, []);

    const clearHover = useCallback(() => {
        setHoveredFeature(null);
        setHoveredVariant(null);
    }, []);

    // ── Render gene arrow polygon ────────────────────────────────────────────
    const renderGeneArrow = useCallback((
        gene: GenomeFeature,
        x1: number,
        x2: number,
        y: number,
    ) => {
        const w = x2 - x1;
        const h = GENE_HEIGHT;
        const yTop = y - h / 2;
        const isForward = gene.strand === 1;
        const clippedLeft = gene.start < viewStart;
        const clippedRight = gene.end > viewEnd;
        const drawHead = isForward ? !clippedRight : !clippedLeft;
        const arrowHead = drawHead ? Math.min(ARROW_SIZE, w) : 0;
        const bodyW = Math.max(0, w - arrowHead);

        let points: string;
        if (isForward) {
            if (drawHead) {
                points = `${x1},${yTop} ${x1 + bodyW},${yTop} ${x2},${y} ${x1 + bodyW},${yTop + h} ${x1},${yTop + h}`;
            } else {
                points = `${x1},${yTop} ${x2},${yTop} ${x2},${yTop + h} ${x1},${yTop + h}`;
            }
        } else {
            if (drawHead) {
                points = `${x1 + arrowHead},${yTop} ${x2},${yTop} ${x2},${yTop + h} ${x1 + arrowHead},${yTop + h} ${x1},${y}`;
            } else {
                points = `${x1},${yTop} ${x2},${yTop} ${x2},${yTop + h} ${x1},${yTop + h}`;
            }
        }

        const fill = isForward ? '#2e7d32' : '#ef6c00';
        const isSelected = selectedVariantPos != null &&
            selectedVariantPos >= gene.start && selectedVariantPos <= gene.end;

        // Show label if gene is wide enough
        const showLabel = w > 30;
        const labelText = gene.name || gene.locus_tag || '';

        return (
            <g
                key={`${gene.locus_tag || gene.name}-${gene.start}`}
                style={{ cursor: 'pointer' }}
                onMouseEnter={(e) => handleFeatureHover(gene, e)}
                onMouseMove={(e) => handleFeatureHover(gene, e)}
                onMouseLeave={clearHover}
                onClick={(e) => {
                    e.stopPropagation();
                    onFeatureClick?.(gene);
                }}
            >
                <polygon
                    points={points}
                    fill={fill}
                    stroke={isSelected ? '#000' : fill}
                    strokeWidth={isSelected ? 2 : 0.5}
                    opacity={0.75}
                />
                {showLabel && (
                    <text
                        x={x1 + w / 2}
                        y={yTop - 3}
                        textAnchor="middle"
                        fontSize={10}
                        fontWeight="bold"
                        fill="#2c3e50"
                        style={{ pointerEvents: 'none', userSelect: 'none' }}
                    >
                        {labelText.length > Math.floor(w / 7) ? labelText.slice(0, Math.floor(w / 7)) + '…' : labelText}
                    </text>
                )}
            </g>
        );
    }, [viewStart, viewEnd, selectedVariantPos, handleFeatureHover, clearHover, onFeatureClick]);

    // ── Ticks ────────────────────────────────────────────────────────────────
    const ticks = useMemo(() => generateTicks(viewStart, viewEnd), [viewStart, viewEnd]);

    // ── Zoom controls ────────────────────────────────────────────────────────
    const panLeft = useCallback(() => {
        const shift = viewSpan * 0.2;
        const newStart = Math.max(0, viewStart - shift);
        const newEnd = newStart + viewSpan;
        onViewChange(Math.round(newStart), Math.round(Math.min(maxPos, newEnd)));
    }, [viewStart, viewSpan, maxPos, onViewChange]);

    const panRight = useCallback(() => {
        const shift = viewSpan * 0.2;
        const newEnd = Math.min(maxPos, viewEnd + shift);
        const newStart = newEnd - viewSpan;
        onViewChange(Math.round(Math.max(0, newStart)), Math.round(newEnd));
    }, [viewEnd, viewSpan, maxPos, onViewChange]);

    const zoomIn = useCallback(() => {
        const center = (viewStart + viewEnd) / 2;
        const newSpan = Math.max(200, viewSpan / 1.5);
        onViewChange(
            Math.round(Math.max(0, center - newSpan / 2)),
            Math.round(Math.min(maxPos, center + newSpan / 2))
        );
    }, [viewStart, viewEnd, viewSpan, maxPos, onViewChange]);

    const zoomOut = useCallback(() => {
        const center = (viewStart + viewEnd) / 2;
        const newSpan = Math.min(maxPos, viewSpan * 1.5);
        onViewChange(
            Math.round(Math.max(0, center - newSpan / 2)),
            Math.round(Math.min(maxPos, center + newSpan / 2))
        );
    }, [viewStart, viewEnd, viewSpan, maxPos, onViewChange]);

    return (
        <div style={{ marginTop: '30px' }}>
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '8px' }}>
                <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
                    <h4 style={{ color: '#2c3e50', margin: 0 }}>Genome Annotation</h4>
                    <span style={{ fontSize: '0.85em', color: '#666' }}>
                        {formatBp(viewStart)} – {formatBp(viewEnd)} ({formatBp(viewSpan)} window)
                    </span>
                </div>
                <div style={{ display: 'flex', gap: '4px' }}>
                    <button onClick={panLeft} title="Pan Left" style={btnStyle}>{'<<'}</button>
                    <button onClick={zoomOut} title="Zoom Out" style={btnStyle}>−</button>
                    <button onClick={zoomIn} title="Zoom In" style={btnStyle}>+</button>
                    <button onClick={panRight} title="Pan Right" style={btnStyle}>{'>>'}</button>
                </div>
            </div>
            <p style={{ fontSize: '0.85em', color: '#888', margin: '0 0 6px 0' }}>
                Drag to pan · Scroll to zoom · Click genes to edit · Click variants for evidence
            </p>

            <div
                ref={containerRef}
                style={{ position: 'relative', width: '100%', userSelect: 'none' }}
            >
                <svg
                    ref={svgRef}
                    width={svgWidth}
                    height={trackHeight + MARGIN.top + MARGIN.bottom}
                    style={{
                        cursor: isDragging ? 'grabbing' : 'grab',
                        display: 'block',
                        border: '1px solid #e0e0e0',
                        borderRadius: '4px',
                        backgroundColor: '#fafafa',
                    }}
                    onWheel={handleWheel}
                    onMouseDown={handleMouseDown}
                    onMouseMove={handleMouseMove}
                    onMouseUp={handleMouseUp}
                >
                    {/* Background */}
                    <rect x={MARGIN.left} y={MARGIN.top} width={plotWidth} height={trackHeight}
                        fill="#fff" stroke="#e0e0e0" strokeWidth={1} />

                    {/* Strand labels */}
                    <text x={MARGIN.left - 8} y={MARGIN.top + FORWARD_Y} textAnchor="end"
                        fontSize={11} fill="#2e7d32" fontWeight="bold" dominantBaseline="middle">+</text>
                    <text x={MARGIN.left - 8} y={MARGIN.top + REVERSE_Y} textAnchor="end"
                        fontSize={11} fill="#ef6c00" fontWeight="bold" dominantBaseline="middle">−</text>
                    <text x={MARGIN.left - 8} y={MARGIN.top + VARIANT_Y} textAnchor="end"
                        fontSize={9} fill="#666" dominantBaseline="middle">Var</text>

                    {/* Center separator line */}
                    <line
                        x1={MARGIN.left} y1={MARGIN.top + VARIANT_Y}
                        x2={MARGIN.left + plotWidth} y2={MARGIN.top + VARIANT_Y}
                        stroke="#e0e0e0" strokeWidth={1} strokeDasharray="4,4"
                    />

                    {/* Clip path for content area */}
                    <defs>
                        <clipPath id="genome-clip">
                            <rect x={MARGIN.left} y={MARGIN.top} width={plotWidth} height={trackHeight} />
                        </clipPath>
                    </defs>

                    <g clipPath="url(#genome-clip)">
                        {/* Offset for margin */}
                        <g transform={`translate(0, ${MARGIN.top})`}>

                            {/* Deletion spans */}
                            {visibleDeletions.map((del, i) => {
                                const x1 = bpToX(Math.max(del.pos, viewStart));
                                const x2 = bpToX(Math.min(del.end || del.pos, viewEnd));
                                return (
                                    <g key={`del-${i}`}>
                                        <rect
                                            x={x1} y={VARIANT_Y - 15} width={Math.max(2, x2 - x1)} height={30}
                                            fill="#ffEBEE" stroke="#ef5350" strokeWidth={1}
                                            strokeDasharray="3,3" opacity={0.8}
                                            style={{ cursor: 'pointer' }}
                                            onClick={(e) => {
                                                e.stopPropagation();
                                                if (del.evidence_links?.length && onEvidenceOpen) {
                                                    onEvidenceOpen(del.evidence_links[0].url);
                                                }
                                            }}
                                            onMouseEnter={(ev) => handleVariantHover(del, ev)}
                                            onMouseMove={(ev) => handleVariantHover(del, ev)}
                                            onMouseLeave={clearHover}
                                        />
                                        <text
                                            x={(x1 + x2) / 2} y={VARIANT_Y - 18}
                                            textAnchor="middle" fontSize={10} fontWeight="bold" fill="#c62828"
                                            style={{ pointerEvents: 'none' }}
                                        >
                                            DEL (−{del.size}bp)
                                        </text>
                                    </g>
                                );
                            })}

                            {/* Variant markers (SNPs, INS) */}
                            {visibleVariants.filter(v => !v.type.includes('DEL')).map((v, i) => {
                                const cx = bpToX(v.pos);
                                const isSelected = v.pos === selectedVariantPos;
                                return (
                                    <g
                                        key={`var-${i}`}
                                        style={{ cursor: 'pointer' }}
                                        onClick={(e) => {
                                            e.stopPropagation();
                                            if (v.evidence_links?.length && onEvidenceOpen) {
                                                onEvidenceOpen(v.evidence_links[0].url);
                                            }
                                            onVariantClick?.(v);
                                        }}
                                        onMouseEnter={(ev) => handleVariantHover(v, ev)}
                                        onMouseMove={(ev) => handleVariantHover(v, ev)}
                                        onMouseLeave={clearHover}
                                    >
                                        {v.type.includes('INS') ? (
                                            // Triangle for insertions
                                            <>
                                                <path
                                                    d={`M${cx},${VARIANT_Y - 6} L${cx + 6},${VARIANT_Y + 6} L${cx - 6},${VARIANT_Y + 6} Z`}
                                                    fill={v.fill}
                                                    stroke={isSelected ? '#000' : 'none'}
                                                    strokeWidth={isSelected ? 2 : 0}
                                                />
                                                {(() => {
                                                    const seqMatch = v.mutation?.match(/\+([ACGTN]+)/i);
                                                    const label = seqMatch && seqMatch[1].length < 6 ? seqMatch[0]
                                                        : v.size > 0 ? `+${v.size}bp` : null;
                                                    return label ? (
                                                        <text x={cx} y={VARIANT_Y - 10} textAnchor="middle"
                                                            fontSize={9} fill="#333" fontWeight="bold"
                                                            style={{ pointerEvents: 'none' }}>{label}</text>
                                                    ) : null;
                                                })()}
                                            </>
                                        ) : (
                                            // Circle for SNP/SUB
                                            <circle
                                                cx={cx} cy={VARIANT_Y} r={5}
                                                fill={v.fill}
                                                stroke={isSelected ? '#000' : 'none'}
                                                strokeWidth={isSelected ? 2 : 0}
                                            />
                                        )}
                                    </g>
                                );
                            })}

                            {/* Forward strand genes */}
                            {forwardGenes.map(({ gene, x1, x2, y }) =>
                                renderGeneArrow(gene, x1, x2, y)
                            )}

                            {/* Reverse strand genes */}
                            {reverseGenes.map(({ gene, x1, x2, y }) =>
                                renderGeneArrow(gene, x1, x2, y)
                            )}

                            {/* Selected variant position indicator */}
                            {selectedVariantPos != null && selectedVariantPos >= viewStart && selectedVariantPos <= viewEnd && (
                                <line
                                    x1={bpToX(selectedVariantPos)} y1={0}
                                    x2={bpToX(selectedVariantPos)} y2={trackHeight}
                                    stroke="#d32f2f" strokeWidth={1.5} strokeDasharray="4,3" opacity={0.5}
                                />
                            )}
                        </g>
                    </g>

                    {/* X Axis */}
                    <g transform={`translate(0, ${MARGIN.top + trackHeight})`}>
                        <line x1={MARGIN.left} y1={0} x2={MARGIN.left + plotWidth} y2={0}
                            stroke="#999" strokeWidth={1} />
                        {ticks.map(t => {
                            const x = bpToX(t);
                            if (x < MARGIN.left || x > MARGIN.left + plotWidth) return null;
                            return (
                                <g key={t}>
                                    <line x1={x} y1={0} x2={x} y2={5} stroke="#999" strokeWidth={1} />
                                    <text x={x} y={18} textAnchor="middle" fontSize={10} fill="#666">
                                        {formatBp(t)}
                                    </text>
                                </g>
                            );
                        })}
                        <text x={MARGIN.left + plotWidth / 2} y={28} textAnchor="middle"
                            fontSize={10} fill="#999">Position (bp)</text>
                    </g>
                </svg>

                {/* Tooltip */}
                {(hoveredFeature || hoveredVariant) && (
                    <div style={{
                        position: 'absolute',
                        left: Math.min(tooltipPos.x + 12, svgWidth - 260),
                        top: tooltipPos.y - 10,
                        backgroundColor: '#fff',
                        border: '1px solid #ccc',
                        borderRadius: '4px',
                        padding: '8px 10px',
                        boxShadow: '0 2px 8px rgba(0,0,0,0.15)',
                        pointerEvents: 'none',
                        zIndex: 10,
                        maxWidth: 250,
                        fontSize: '0.85em',
                    }}>
                        {hoveredFeature && (
                            <>
                                <div style={{ fontWeight: 'bold', color: '#2e7d32' }}>
                                    {hoveredFeature.name || hoveredFeature.locus_tag}
                                </div>
                                <div style={{ color: '#666', marginTop: 2 }}>{hoveredFeature.product}</div>
                                <div style={{ color: '#999', marginTop: 2, fontSize: '0.9em' }}>
                                    {hoveredFeature.start.toLocaleString()} – {hoveredFeature.end.toLocaleString()} ({hoveredFeature.strand === 1 ? '+' : '−'})
                                </div>
                            </>
                        )}
                        {hoveredVariant && (
                            <>
                                <div style={{ fontWeight: 'bold' }}>
                                    {hoveredVariant.type} @ {hoveredVariant.pos.toLocaleString()}
                                </div>
                                <div style={{ color: '#d32f2f', marginTop: 2 }}>
                                    {hoveredVariant.mutation} ({hoveredVariant.freq}%)
                                </div>
                                {hoveredVariant.gene && (
                                    <div style={{ fontStyle: 'italic', marginTop: 2 }}>{hoveredVariant.gene}</div>
                                )}
                                {hoveredVariant.annotation && (
                                    <div style={{ color: '#666', marginTop: 2, fontSize: '0.9em' }}>{hoveredVariant.annotation}</div>
                                )}
                                {hoveredVariant.size > 1 && (
                                    <div style={{ marginTop: 2, fontWeight: 'bold', fontSize: '0.9em' }}>Size: {hoveredVariant.size} bp</div>
                                )}
                            </>
                        )}
                    </div>
                )}
            </div>
        </div>
    );
};

const btnStyle: React.CSSProperties = {
    padding: '5px 10px',
    cursor: 'pointer',
    fontWeight: 'bold',
    border: '1px solid #ccc',
    borderRadius: '4px',
    backgroundColor: '#f5f5f5',
    fontSize: '0.9em',
};

export default GenomeBrowser;
