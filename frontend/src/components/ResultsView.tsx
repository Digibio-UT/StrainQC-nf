import React, { useEffect, useState, useMemo } from 'react';
import axios from 'axios';
import { 
    ComposedChart, 
    Area, 
    Scatter, 
    ScatterChart,
    XAxis, 
    YAxis, 
    CartesianGrid, 
    Tooltip, 
    Legend, 
    ResponsiveContainer,
    Cell,
    ReferenceArea
} from 'recharts';
import GenomeBrowser from './GenomeBrowser';
import CopyNumberPanel from './CopyNumberPanel';
import type { ChromosomeCopyNumber } from './CopyNumberPanel';

// Removed unused interface AnalysisStatus


// Mapping for Yeast Chromosome Accessions
const ACCESSION_MAP: Record<string, string> = {
    "BK006935": "Chr I",
    "BK006936": "Chr II",
    "BK006937": "Chr III",
    "BK006938": "Chr IV",
    "BK006939": "Chr V",
    "BK006940": "Chr VI",
    "BK006941": "Chr VII",
    "BK006934": "Chr VIII",
    "BK006942": "Chr IX",
    "BK006943": "Chr X",
    "BK006944": "Chr XI",
    "BK006945": "Chr XII",
    "BK006946": "Chr XIII",
    "BK006947": "Chr XIV",
    "BK006948": "Chr XV",
    "BK006949": "Chr XVI",
    "AJ011856": "Mitochondria",
    "NC_001133": "Chr I", "NC_001134": "Chr II", "NC_001135": "Chr III",
    "NC_001136": "Chr IV", "NC_001137": "Chr V", "NC_001138": "Chr VI",
    "NC_001139": "Chr VII", "NC_001140": "Chr VIII", "NC_001141": "Chr IX",
    "NC_001142": "Chr X", "NC_001143": "Chr XI", "NC_001144": "Chr XII",
    "NC_001145": "Chr XIII", "NC_001146": "Chr XIV", "NC_001147": "Chr XV",
    "NC_001148": "Chr XVI", "NC_001224": "Mitochondria"
};

const getDisplayName = (id: string) => {
    // Strip version if present (e.g. .1, .2)
    const baseId = id.split('.')[0];
    return ACCESSION_MAP[baseId] || id;
};

interface ResultsViewProps {
    jobId: string;
}

interface Mutation {
    type: string;
    seq_id: string;
    position: string;
    mutation: string;
    freq: string;
    annotation: string;
    gene: string;
    description: string;
    details: string[];
    evidence_links?: { text: string; url: string }[];
}

interface Contaminant {
    name: string;
    rank: string;
    percentage: number;
    reads: number;
}

interface BlastHit {
    title: string;
    count: number;
    identity: number;
    evalue: number;
    queries: string;
}

interface MashHit {
    identity: number;
    shared_hashes: string;
    median_multiplicity: number;
    p_value: number;
    query_id: string;
    comment: string;
}

interface AnalysisResult {
    job_id: string;
    status: string;
    variant_caller?: string;
    polymorphism_mode?: boolean;
    mutations?: Mutation[];
    contaminants?: Contaminant[];
    blast_results?: BlastHit[];
    mash_results?: MashHit[];
    mash_unmapped_results?: MashHit[];
    snp_blast_results?: BlastHit[];
    snp_query_file?: string;
    hc_snps_file?: string;
    coverage_file?: string;
    unmapped_files?: string[];
    // CNVkit copy number (gatk path, only when run with --parent_sample)
    chromosome_copy_number?: ChromosomeCopyNumber[];
    aneuploidy?: { seq_id: string; call: string; mean_cn: number }[];
    cnv_tool?: string;
    read_stats?: {
        total_reads: number;
        mapped_percentage: number;
        unmapped_percentage: number;
        unmapped_reads: number;
    };
}

interface CoverageData {
    chrom: string;
    pos: number;
    depth: number;
}

interface GenomeFeature {
    start: number;
    end: number;
    strand: number;
    name: string;
    locus_tag: string;
    product: string;
    seq_id?: string;
    parts?: { start: number; end: number }[];
}

// --- Custom Shape Components defined outside to satisfy Recharts ---

// GeneArrowShape removed — replaced by SVG-based GenomeBrowser component

// Custom shape for rendering variants (Insertions/SNPs)
const VariantShape = (props: any) => {
    const { cx, cy, payload } = props;
    // Check for selected variant via payload.pos or strictly strictly if passed directly
    // Note: Scatter rendering passes payload which allows checking pos
    const isSelected = payload.pos === payload.selectedVariantPos; 
    const stroke = isSelected ? '#000' : 'none';
    const strokeWidth = isSelected ? 2 : 0;
    
    if (payload.type.includes('INS')) {
        // Draw a triangle pointing up for Insertions
        let label = null;
        // Check for sequence insertion like +G, +ACGT with strict matching or extracting
        const seqMatch = payload.mutation ? payload.mutation.match(/\+([ACGTN]+)/i) : null;
        
        if (seqMatch && seqMatch[1].length < 6) {
             // If it's a short sequence, show the sequence itself e.g., +G
             label = seqMatch[0];
        } else if (payload.size > 0) {
             // Otherwise show the size e.g., +150bp
             // Note: payload.size includes length of sequence if it was parsed
             label = `+${payload.size}bp`;
        }

        return (
            <g>
                <path d={`M${cx},${cy-6} L${cx+6},${cy+6} L${cx-6},${cy+6} Z`} fill={payload.fill} stroke={stroke} strokeWidth={strokeWidth} />
                {label && (
                    <text x={cx} y={cy-10} textAnchor="middle" fontSize={10} fill="#333" fontWeight="bold">
                        {label}
                    </text>
                )}
            </g>
        );
    }

    if (payload.type.includes('DEL')) {
        // Draw an X mark for deletions
        const s = 5; // half-size
        let label = null;
        const deltaMatch = payload.mutation ? payload.mutation.match(/[Δ∆-]\s*([\d,]+)/) : null;
        if (deltaMatch) {
            label = `Δ${deltaMatch[1]}bp`;
        }
        return (
            <g>
                <line x1={cx-s} y1={cy-s} x2={cx+s} y2={cy+s} stroke={payload.fill} strokeWidth={2} />
                <line x1={cx+s} y1={cy-s} x2={cx-s} y2={cy+s} stroke={payload.fill} strokeWidth={2} />
                {label && (
                    <text x={cx} y={cy-10} textAnchor="middle" fontSize={10} fill="#333" fontWeight="bold">
                        {label}
                    </text>
                )}
            </g>
        );
    }
    
    // Default circle for SNP/SUB
    return <circle cx={cx} cy={cy} r={5} fill={payload.fill} stroke={stroke} strokeWidth={strokeWidth} />;
};


// Helper to sort chromosomes: Chr I, Chr II ... Mitochondria (last)
// Assumes getDisplayName returns simplified names like "Chr I", "Mitochondria"
const sortChromosomes = (chroms: string[]) => {
    // Custom sort order logic
    
    // We can rely on getDisplayName map.
    // However, getDisplayName uses the ID.
    // Let's implement a comparator.
    
    // Simple helper to convert Roman to integer for proper sorting
    const romanToInt = (s: string): number => {
        const romanMap: {[key: string]: number} = { 
            'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 
            'VI': 6, 'VII': 7, 'VIII': 8, 'IX': 9, 'X': 10, 
            'XI': 11, 'XII': 12, 'XIII': 13, 'XIV': 14, 'XV': 15, 'XVI': 16 
        };
        return romanMap[s] || 999;
    };

    return [...chroms].sort((a, b) => { // Create a copy before sorting
        const dispA = getDisplayName(a);
        const dispB = getDisplayName(b);
        
        // Ensure Mitochondria is last
        if (dispA.includes("Mitochondria")) return 1;
        if (dispB.includes("Mitochondria")) return -1;
        
        // Try to parse Roman numerals or Numbers from "Chr X"
        // Regex for "Chr <Roman/Number>"
        const matchA = dispA.match(/Chr\s+([IVXLCDM]+|\d+)/);
        const matchB = dispB.match(/Chr\s+([IVXLCDM]+|\d+)/);
        
        let valA: number | string = dispA;
        let valB: number | string = dispB;

        if (matchA) {
            if (/^\d+$/.test(matchA[1])) valA = parseInt(matchA[1]);
            else valA = romanToInt(matchA[1]);
        }
        
        if (matchB) {
            if (/^\d+$/.test(matchB[1])) valB = parseInt(matchB[1]);
            else valB = romanToInt(matchB[1]);
        }
        
        if (typeof valA === 'number' && typeof valB === 'number') return valA - valB;
        if (typeof valA === 'number') return -1; // Numbers first
        if (typeof valB === 'number') return 1;
        
        // Fallback to string comparison
        return dispA.localeCompare(dispB);
    });
};

const ResultsView: React.FC<ResultsViewProps> = ({ jobId }) => {
    const [status, setStatus] = useState<string>('loading');
    const [processingSteps, setProcessingSteps] = useState<Record<string, string>>({});
    const [results, setResults] = useState<AnalysisResult | null>(null);
    const [coverage, setCoverage] = useState<CoverageData[]>([]);
    const [features, setFeatures] = useState<GenomeFeature[]>([]);
    const [error, setError] = useState<string>('');
    const [selectedVariantPos, setSelectedVariantPos] = useState<number | null>(null);
    const [chromosomes, setChromosomes] = useState<string[]>([]);
    const [selectedChrom, setSelectedChrom] = useState<string>('');
    
    // Top Plot State (Overview)
    const [zoomWindowTop, setZoomWindowTop] = useState<number>(5000); 
    const [viewCenterTop, setViewCenterTop] = useState<number | null>(null);

    // Bottom Plot State (Detailed)
    const [zoomWindowBottom, setZoomWindowBottom] = useState<number>(5000); 
    const [viewCenterBottom, setViewCenterBottom] = useState<number | null>(null);

    const [evidenceUrl, setEvidenceUrl] = useState<string | null>(null);
    const [consensusThreshold, setConsensusThreshold] = useState<string>("0.9");
    const [minQualThreshold, setMinQualThreshold] = useState<string>("0");

    // --- Consensus Preview State ---
    const [showConsensusModal, setShowConsensusModal] = useState(false);
    const [consensusProposals, setConsensusProposals] = useState<any[]>([]);
    
    // --- Editing State ---
    const [editingFeature, setEditingFeature] = useState<GenomeFeature | null>(null);
    const [isAddModalOpen, setIsAddModalOpen] = useState(false);
    const [editForm, setEditForm] = useState({ name: '', locus_tag: '', product: '' });
    const [addForm, setAddForm] = useState({ name: '', locus_tag: '', product: '', start: '', end: '', strand: '1', seq_id: '' });

    // --- Gene Search State ---
    const [searchQuery, setSearchQuery] = useState('');
    const [searchResults, setSearchResults] = useState<GenomeFeature[]>([]);
    const [showSearchResults, setShowSearchResults] = useState(false);
    const searchRef = React.useRef<HTMLDivElement>(null);

    // Close search dropdown on outside click
    React.useEffect(() => {
        const handleClickOutside = (e: MouseEvent) => {
            if (searchRef.current && !searchRef.current.contains(e.target as Node)) {
                setShowSearchResults(false);
            }
        };
        document.addEventListener('mousedown', handleClickOutside);
        return () => document.removeEventListener('mousedown', handleClickOutside);
    }, []);

    const handleSearch = (query: string) => {
        setSearchQuery(query);
        if (query.length < 2) {
            setSearchResults([]);
            setShowSearchResults(false);
            return;
        }
        const q = query.toLowerCase();
        const matches = features.filter(f => 
            (f.name && f.name.toLowerCase().includes(q)) ||
            (f.locus_tag && f.locus_tag.toLowerCase().includes(q)) ||
            (f.product && f.product.toLowerCase().includes(q))
        ).slice(0, 20);
        setSearchResults(matches);
        setShowSearchResults(matches.length > 0);
    };

    const handleSelectSearchResult = (feature: GenomeFeature) => {
        // Switch chromosome if needed — reset top plot to full overview
        if (feature.seq_id && feature.seq_id !== selectedChrom) {
            setSelectedChrom(feature.seq_id);
            setViewCenterTop(null);
        }
        // Navigate genome browser (bottom plot) to feature center
        const center = Math.round((feature.start + feature.end) / 2);
        const featureSize = feature.end - feature.start;
        const viewWindow = Math.max(featureSize * 3, 5000);
        setViewCenterBottom(center);
        setZoomWindowBottom(viewWindow);
        setSelectedVariantPos(null);
        // Close search
        setShowSearchResults(false);
        setSearchQuery('');
    };

    // --- Edit Handlers ---
    const openEditModal = (feature: GenomeFeature) => {
        setEditingFeature(feature);
        setEditForm({
            name: feature.name || '',
            locus_tag: feature.locus_tag || '',
            product: feature.product || ''
        });
    };

    const closeEditModal = () => {
        setEditingFeature(null);
        setEditForm({ name: '', locus_tag: '', product: '' });
    };

    const handleSaveEdit = async () => {
        if (!editingFeature) return;
        try {
            // Determine ID (use persistent original_id if available)
            // @ts-ignore
            const originalId = editingFeature.original_id || editingFeature.locus_tag || editingFeature.name;
            await axios.post(`/api/v1/analysis/annotations/${jobId}`, {
                original_id: originalId,
                action: 'edit',
                new_name: editForm.name,
                new_locus_tag: editForm.locus_tag,
                new_product: editForm.product
            });
            closeEditModal();
            // Refresh features
            const res = await axios.get(`/api/v1/analysis/genome_features/${jobId}`);
            setFeatures(res.data.features);
        } catch (e) {
            console.error("Failed to save edit", e);
            alert("Failed to save edit.");
        }
    };

    const handleDeleteFeature = async () => {
        if (!editingFeature) return;
        if (!window.confirm(`Are you sure you want to delete ${editingFeature.name}?`)) return;
        try {
            // @ts-ignore
            const originalId = editingFeature.original_id || editingFeature.locus_tag || editingFeature.name;
            await axios.post(`/api/v1/analysis/annotations/${jobId}`, {
                original_id: originalId,
                action: 'delete'
            });
            closeEditModal();
            // Refresh features
            const res = await axios.get(`/api/v1/analysis/genome_features/${jobId}`);
            setFeatures(res.data.features);
        } catch (e) {
            console.error("Failed to delete feature", e);
            alert("Failed to delete feature.");
        }
    };

    const handleAddFeature = async () => {
        try {
            await axios.post(`/api/v1/analysis/annotations/${jobId}`, {
                original_id: `new_${Date.now()}`, // Temporary ID
                action: 'add',
                new_name: addForm.name,
                new_locus_tag: addForm.locus_tag,
                new_product: addForm.product,
                start: parseInt(addForm.start),
                end: parseInt(addForm.end),
                strand: parseInt(addForm.strand),
                seq_id: addForm.seq_id || chromosomes[0] || 'unknown'
            });
            setIsAddModalOpen(false);
            setAddForm({ name: '', locus_tag: '', product: '', start: '', end: '', strand: '1', seq_id: '' });
            // Refresh features
            const res = await axios.get(`/api/v1/analysis/genome_features/${jobId}`);
            setFeatures(res.data.features);
        } catch (e) {
            console.error("Failed to add feature", e);
            alert("Failed to add feature. Ensure fields are valid.");
        }
    };

    // Process mutations for visualization - optimized
    const variantData = useMemo(() => {
        if (!results?.mutations) return [];
        return results.mutations
            .filter(m => !selectedChrom || m.seq_id === selectedChrom || !m.seq_id)
            .map(m => {
            try {
                const pos = parseInt(m.position.replace(/,/g, ''), 10);
                const freq = parseFloat(m.freq.replace('%', ''));
                
                // Determine color based on mutation type
                let color = '#d32f2f'; // Default red for SNP
                if (m.type.includes('DEL')) color = '#1976d2'; // Blue for Deletions
                if (m.type.includes('INS')) color = '#388e3c'; // Green for Insertions
                if (m.type.includes('SUB')) color = '#fbc02d'; // Yellow for Substitutions
                
                // Parse size/end for Deletions and Insertions
                let size = 1;
                let end = pos;
                
                if (m.type.includes('DEL')) {
                    // Robust parsing for deletion size
                    // Matches: Δ123, Δ 123, 123bp, 123 bp, -123
                    let sizeFound = false;
                    
                    // Pattern 1: Delta or minus followed by number
                    const deltaMatch = m.mutation.match(/[Δ∆-]\s*([\d,]+)/);
                    if (deltaMatch) {
                        size = parseInt(deltaMatch[1].replace(/,/g, ''), 10);
                        sizeFound = true;
                    } 
                    
                    // Pattern 2: Number followed by bp (if not found yet)
                    if (!sizeFound) {
                        const bpMatch = m.mutation.match(/([\d,]+)\s*bp/i);
                        if (bpMatch) {
                            size = parseInt(bpMatch[1].replace(/,/g, ''), 10);
                            sizeFound = true;
                        }
                    }

                    if (sizeFound && size > 0) {
                        end = pos + size;
                    }
                } else if (m.type.includes('INS')) {
                     // Check for "+1,199 bp" format
                    const sizeMatch = m.mutation.match(/\+([\d,]+)\s*bp/);
                    if (sizeMatch) {
                        size = parseInt(sizeMatch[1].replace(/,/g, ''), 10);
                    } else {
                        // Check for sequence "+ACGT"
                        const seqMatch = m.mutation.match(/\+([ACGTN]+)/i);
                        if (seqMatch) {
                            size = seqMatch[1].length;
                        }
                    }
                }

                return {
                    x: pos, 
                    y: freq, 
                    pos,
                    end,
                    size,
                    freq,
                    type: m.type,
                    mutation: m.mutation,
                    gene: m.gene,
                    annotation: m.annotation,
                    evidence_links: m.evidence_links,
                    fill: color // Add color property
                };
            } catch (e) {
                return null;
            }
        }).filter(item => item !== null).sort((a, b) => (a?.pos || 0) - (b?.pos || 0));
    }, [results?.mutations, selectedChrom]);

    const maxPos = useMemo(() => {
        const filteredCov = coverage.filter(c => !selectedChrom || c.chrom === selectedChrom || !c.chrom);
        if (filteredCov.length > 0) {
            // Find the true maximum position in the filtered dataset
            // necessary because if multiple chromosomes are combined, the last point isn't guaranteed to be the max
            return Math.max(...filteredCov.map(c => c.pos));
        }
        return 0;
    }, [coverage, selectedChrom]);

    // Set initial view centers when data loads or chromosome changes
    useEffect(() => {
        if (maxPos > 0) {
            // Top Plot: Full Overview
            if (viewCenterTop === null) {
                setZoomWindowTop(maxPos);
                setViewCenterTop(maxPos / 2);
            }

            // Bottom Plot: Focused on first variant (if any)
            if (viewCenterBottom === null) {
                if (variantData.length > 0) {
                    setViewCenterBottom(variantData[0].pos);
                    setZoomWindowBottom(5000); // 5kb window
                    // Also select the variant without triggering zoom logic (handled separately)
                    if (selectedVariantPos === null) {
                        setSelectedVariantPos(variantData[0].pos);
                    }
                } else {
                    setViewCenterBottom(maxPos / 2);
                    setZoomWindowBottom(5000);
                }
            }
        }
    }, [variantData, selectedVariantPos, maxPos, viewCenterTop, viewCenterBottom]);

    // Handle jump to variant on bottom plot when selectedVariantPos changes
    // This happens when clicking a point on the top plot
    useEffect(() => {
        if (selectedVariantPos !== null) {
            setViewCenterBottom(selectedVariantPos);
            // Maintain current zoom level on bottom plot (do not reset)
        }
    }, [selectedVariantPos]);

    // Prepare coverage data for chart: filter by chromosome and ensure start at 0
    const chartCoverageData = useMemo(() => {
        let data = coverage.filter(c => !selectedChrom || c.chrom === selectedChrom || !c.chrom);
        
        // Safety: ensure data is sorted by position (assuming single chrom or ordered)
        // If multiple chromosomes are shown (All), sorting by pos mixes them if they overlap in range. 
        // But "All" usually implies concatenation or overlay. Overlay is messy.
        // Assuming "All" view implies showing first chrom? Or just raw list?
        // If sorting is needed: data.sort((a,b) => a.pos - b.pos);
        
        // Inject (0,0) point to fix "rectangle artifact" at start
        if (data.length > 0 && data[0].pos > 0) {
            // Create a pseudo start point
            const startPoint = { ...data[0], pos: 0, depth: 0 };
            data = [startPoint, ...data];
        }
        return data;
    }, [coverage, selectedChrom]);



    // Custom shape definitions moved outside component


    // Process features for visualization
    const featureData = useMemo(() => {
        return features
            .filter(f => !selectedChrom || f.seq_id === selectedChrom || !f.seq_id)
            .map(f => ({
            x: (f.start + f.end) / 2, // Center of gene
            y: f.strand, // 1 or -1
            ...f
        }));
    }, [features, selectedChrom]);

    const zoomDomainTop = useMemo(() => {
        if (viewCenterTop === null) return [0, maxPos];
        const start = Math.max(0, viewCenterTop - zoomWindowTop / 2);
        const end = viewCenterTop + zoomWindowTop / 2;
        return [start, end];
    }, [viewCenterTop, zoomWindowTop, maxPos]);

    const zoomDomainBottom = useMemo(() => {
        if (viewCenterBottom === null) return [0, maxPos];
        const start = Math.max(0, viewCenterBottom - zoomWindowBottom / 2);
        const end = viewCenterBottom + zoomWindowBottom / 2;
        return [start, end];
    }, [viewCenterBottom, zoomWindowBottom, maxPos]);

    const CustomTooltip = ({ active, payload, label }: any) => {
        if (active && payload && payload.length) {
            // Check if we are hovering over a scatter point (Variant)
            const variantPayload = payload.find((p: any) => p.name === 'Variants');
            const deletionPayload = payload.find((p: any) => p.name === 'Deletions');
            const coveragePayload = payload.find((p: any) => p.name === 'Coverage');
            const featurePayload = payload.find((p: any) => p.name === 'Genes');

            if (variantPayload || deletionPayload) {
                const data = variantPayload ? variantPayload.payload : deletionPayload.payload;
                return (
                    <div style={{ backgroundColor: '#fff', padding: '10px', border: '1px solid #ccc', borderRadius: '4px' }}>
                        <p style={{ fontWeight: 'bold', margin: 0 }}>{data.type} @ {data.pos}</p>
                        <p style={{ margin: '4px 0', color: '#d32f2f' }}>{data.mutation} ({data.freq}%)</p>
                        <p style={{ margin: 0, fontSize: '0.9em', fontStyle: 'italic' }}>{data.gene}</p>
                        <p style={{ margin: 0, fontSize: '0.8em', color: '#666' }}>{data.annotation}</p>
                        {data.size > 1 && <p style={{ margin: '4px 0 0 0', fontSize: '0.8em', fontWeight: 'bold' }}>Size: {data.size} bp</p>}
                    </div>
                );
            }

            if (featurePayload) {
                const data = featurePayload.payload;
                return (
                    <div style={{ backgroundColor: '#fff', padding: '10px', border: '1px solid #ccc', borderRadius: '4px', maxWidth: '300px' }}>
                        <p style={{ fontWeight: 'bold', margin: 0, color: '#2e7d32' }}>{data.name || data.locus_tag}</p>
                        <p style={{ margin: '4px 0', fontSize: '0.9em' }}>{data.product}</p>
                        <p style={{ margin: 0, fontSize: '0.8em', color: '#666' }}>
                            {data.start.toLocaleString()} - {data.end.toLocaleString()} ({data.strand === 1 ? '+' : '-'})
                        </p>
                    </div>
                );
            }

            if (coveragePayload) {
                return (
                    <div style={{ backgroundColor: '#fff', padding: '10px', border: '1px solid #ccc', borderRadius: '4px' }}>
                        <p style={{ margin: 0 }}>Position: {label}</p>
                        <p style={{ margin: 0, color: '#8884d8' }}>Depth: {coveragePayload.value}</p>
                    </div>
                );
            }
        }
        return null;
    };

    useEffect(() => {
        let intervalId: any;
        setStatus('loading');

        const checkStatus = async () => {
            try {
                const response = await axios.get(`/api/v1/analysis/status/${jobId}`);
                const currentStatus = response.data.status;
                setStatus(currentStatus);
                
                // Update processing steps if available
                if (response.data.steps) {
                    setProcessingSteps(response.data.steps);
                }

                if (currentStatus === 'completed') {
                    clearInterval(intervalId);
                    fetchResults();
                } else if (currentStatus === 'failed') {
                    clearInterval(intervalId);
                    const errorMsg = response.data.error 
                        ? `Analysis failed: ${response.data.error}` 
                        : 'Analysis failed.';
                    setError(errorMsg);
                }
            } catch (err) {
                console.error(err);
                // Don't stop polling immediately on transient errors, but maybe limit retries in real app
            }
        };

        const fetchResults = async () => {
            try {
                const response = await axios.get(`/api/v1/analysis/results/${jobId}`);
                setResults(response.data);
                
                // Extract unique Chromosomes
                if (response.data.mutations) {
                    const muts: Mutation[] = response.data.mutations;
                    const chroms = sortChromosomes(Array.from(new Set(
                        muts.map(m => m.seq_id).filter(Boolean)
                    )));
                    
                    if (chroms.length > 0) {
                        setChromosomes(chroms);
                        setSelectedChrom(chroms[0]);
                    }
                }
                
                if (response.data.coverage_file) {
                    fetchCoverage();
                    fetchFeatures();
                }
            } catch (err) {
                setError('Failed to fetch results.');
                console.error(err);
            }
        };
        
        const fetchCoverage = async (chrom = "") => {
            try {
                // If fetching a specific chromosome, increase resolution or check cache
                // But for now, just fetch with param
                const url = chrom 
                    ? `/api/v1/analysis/coverage/${jobId}?resolution=2000&chrom=${chrom}`
                    : `/api/v1/analysis/coverage/${jobId}?resolution=2000`;
                
                const response = await axios.get(url);
                const coverageData: CoverageData[] = response.data.data;
                setCoverage(coverageData);

                // Extract chromosomes from coverage
                if (!chrom && coverageData.length > 0) {
                     const covChroms = Array.from(new Set(coverageData.map((c: any) => c.chrom).filter(Boolean)));
                     
                     setChromosomes(prev => {
                         const combined = sortChromosomes(Array.from(new Set([...prev, ...covChroms])));
                         // If no chromosome was selected yet (e.g. no mutations found), select the first one from coverage
                         if (prev.length === 0 && combined.length > 0) {
                             setTimeout(() => setSelectedChrom(combined[0]), 0);
                         }
                         return combined;
                     });
                }
            } catch (err) {
                console.error("Failed to fetch coverage data", err);
            }
        };

        const fetchFeatures = async () => {
            try {
                const response = await axios.get(`/api/v1/analysis/genome_features/${jobId}`);
                const feats = response.data.features;
                setFeatures(feats);

                // Populate chromosomes from features if not already set (e.g. if coverage failed but features exist)
                if (feats.length > 0) {
                     const featChroms = Array.from(new Set(feats.map((f: any) => f.seq_id))).filter(Boolean) as string[];
                     
                     if (featChroms.length > 0) {
                        setChromosomes(prev => {
                            const unique = new Set([...prev, ...featChroms]);
                            const combined = sortChromosomes(Array.from(unique));
                            // Select first chromosome if none selected
                            if (prev.length === 0 && combined.length > 0) {
                                setTimeout(() => setSelectedChrom(combined[0]), 0);
                            }
                            return combined;
                        });
                     }
                }
            } catch (err) {
                console.error("Failed to fetch genome features", err);
            }
        };

        // Initial check
        checkStatus();
        
        // Poll every 2 seconds
        intervalId = setInterval(checkStatus, 2000);

        return () => clearInterval(intervalId);
    }, [jobId]);

    // Refetch coverage when chromosome changes
    // But we need to be careful not to cause loops or wipe out initial full data if needed for dropdown
    useEffect(() => {
        if (selectedChrom && status === 'completed') {
            const fetchChromCoverage = async () => {
                try {
                     const response = await axios.get(`/api/v1/analysis/coverage/${jobId}?resolution=2000&chrom=${selectedChrom}`);
                     setCoverage(response.data.data);
                } catch (e) { console.error(e); }
            };
            fetchChromCoverage();
        }
    }, [selectedChrom, jobId, status]);

    if (error) {
        return <div style={{ color: 'red', marginTop: '20px' }}>{error}</div>;
    }

    if (status !== 'completed') {
        const stepItems = Object.entries(processingSteps).map(([step, state]) => (
            <li key={step} style={{ textAlign: 'left', marginLeft: '30%' }}>
                <strong>{step}:</strong> {state}
            </li>
        ));

        return (
            <div style={{ marginTop: '20px', padding: '20px', backgroundColor: '#eef', borderRadius: '8px' }}>
                <h3>Job Running</h3>
                <p style={{ fontWeight: 'bold' }}>Current Status: {status.toUpperCase()}</p>
                
                {stepItems.length > 0 && (
                    <div style={{ marginTop: '15px', marginBottom: '15px' }}>
                        <p>Processing Steps:</p>
                        <ul style={{ listStyleType: 'none', padding: 0 }}>
                            {stepItems}
                        </ul>
                    </div>
                )}
                
                <div className="loader" style={{ border: '4px solid #f3f3f3', borderTop: '4px solid #3498db', borderRadius: '50%', width: '30px', height: '30px', animation: 'spin 2s linear infinite', margin: '0 auto' }}></div>
                <style>{`@keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }`}</style>
            </div>
        );
    }

    return (
        <div style={{ marginTop: '30px', textAlign: 'left' }}>
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
                    <h2>Analysis Results</h2>
                    {results?.variant_caller && (
                        <span style={{
                            fontSize: '0.75em',
                            padding: '3px 8px',
                            borderRadius: '12px',
                            backgroundColor: results.variant_caller === 'gatk' ? '#e8f5e9' : '#e3f2fd',
                            color: results.variant_caller === 'gatk' ? '#2e7d32' : '#1565c0',
                            border: `1px solid ${results.variant_caller === 'gatk' ? '#a5d6a7' : '#90caf9'}`,
                            fontWeight: 600,
                        }}>
                            {results.variant_caller === 'gatk' 
                                ? (results.polymorphism_mode ? 'GATK Mutect2' : 'GATK') 
                                : 'breseq'}
                        </span>
                    )}
                </div>
                <div style={{ display: 'flex', alignItems: 'center', gap: '15px' }}>
                    
                    {/* Chromosome Selector */}
                    {chromosomes.length > 0 && (
                        <div style={{ display: 'flex', alignItems: 'center', gap: '8px', backgroundColor: '#e3f2fd', padding: '5px 10px', borderRadius: '5px', border: '1px solid #90caf9' }}>
                            <label style={{ fontSize: '0.85em', color: '#1565c0', fontWeight: 'bold' }}>
                                Chromosome:
                            </label>
                            <select 
                                value={selectedChrom}
                                onChange={(e) => {
                                    setSelectedChrom(e.target.value);
                                    setSelectedVariantPos(null);
                                    setViewCenterTop(null);
                                    setViewCenterBottom(null);
                                    setZoomWindowTop(5000);
                                    setZoomWindowBottom(5000);
                                }}
                                style={{ padding: '4px', borderRadius: '4px', border: '1px solid #ccc' }}
                            >
                                <option value="" disabled>Select Chromosome</option>
                                {chromosomes.map(chrom => (
                                    <option key={chrom} value={chrom}>{getDisplayName(chrom)} ({chrom})</option>
                                ))}
                            </select>
                        </div>
                    )}

                    {/* Gene Search */}
                    {features.length > 0 && (
                        <div ref={searchRef} style={{ position: 'relative', display: 'flex', alignItems: 'center', gap: '6px', backgroundColor: '#fff8e1', padding: '5px 10px', borderRadius: '5px', border: '1px solid #ffe082' }}>
                            <label style={{ fontSize: '0.85em', color: '#f57f17', fontWeight: 'bold', whiteSpace: 'nowrap' }}>
                                🔍
                            </label>
                            <input
                                type="text"
                                placeholder="Search gene, locus, product..."
                                value={searchQuery}
                                onChange={(e) => handleSearch(e.target.value)}
                                onFocus={() => { if (searchResults.length > 0) setShowSearchResults(true); }}
                                style={{ padding: '4px 8px', borderRadius: '4px', border: '1px solid #ccc', width: '200px', fontSize: '0.85em' }}
                            />
                            {showSearchResults && (
                                <div style={{
                                    position: 'absolute',
                                    top: '100%',
                                    left: 0,
                                    right: 0,
                                    minWidth: '350px',
                                    maxHeight: '300px',
                                    overflowY: 'auto',
                                    backgroundColor: '#fff',
                                    border: '1px solid #ccc',
                                    borderRadius: '4px',
                                    boxShadow: '0 4px 12px rgba(0,0,0,0.15)',
                                    zIndex: 1000,
                                }}>
                                    {searchResults.map((f, i) => (
                                        <div
                                            key={`${f.seq_id}_${f.start}_${i}`}
                                            onClick={() => handleSelectSearchResult(f)}
                                            style={{
                                                padding: '8px 12px',
                                                cursor: 'pointer',
                                                borderBottom: i < searchResults.length - 1 ? '1px solid #eee' : 'none',
                                                fontSize: '0.85em',
                                            }}
                                            onMouseEnter={(e) => (e.currentTarget.style.backgroundColor = '#f5f5f5')}
                                            onMouseLeave={(e) => (e.currentTarget.style.backgroundColor = '#fff')}
                                        >
                                            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'baseline' }}>
                                                <span style={{ fontWeight: 600 }}>{f.name || f.locus_tag || 'unnamed'}</span>
                                                <span style={{ fontSize: '0.8em', color: '#888', marginLeft: '8px' }}>
                                                    {getDisplayName(f.seq_id || '')} : {f.start.toLocaleString()}–{f.end.toLocaleString()}
                                                </span>
                                            </div>
                                            {f.product && (
                                                <div style={{ fontSize: '0.8em', color: '#666', marginTop: '2px' }}>{f.product}</div>
                                            )}
                                            {f.locus_tag && f.name && (
                                                <div style={{ fontSize: '0.75em', color: '#999', marginTop: '1px' }}>{f.locus_tag}</div>
                                            )}
                                        </div>
                                    ))}
                                </div>
                            )}
                        </div>
                    )}

                    <div style={{ display: 'flex', alignItems: 'center', gap: '8px', backgroundColor: '#f5f5f5', padding: '5px 10px', borderRadius: '5px' }}>
                        <label style={{ fontSize: '0.85em', color: '#666', fontWeight: 500 }}>
                            Mutation Threshold:
                        </label>
                        <input 
                            type="number" 
                            step="0.05" 
                            min="0.0" 
                            max="1.0" 
                            value={consensusThreshold} 
                            onChange={(e) => setConsensusThreshold(e.target.value)}
                            title="Minimum frequency (0.0-1.0) to include a mutation in the consensus genome"
                            style={{ 
                                padding: '6px', 
                                borderRadius: '4px', 
                                border: '1px solid #ccc', 
                                width: '60px',
                                textAlign: 'center'
                            }} 
                        />
                    </div>
                    {results?.variant_caller === 'gatk' && (
                        <div style={{ display: 'flex', alignItems: 'center', gap: '8px', backgroundColor: '#f5f5f5', padding: '5px 10px', borderRadius: '5px' }}>
                            <label style={{ fontSize: '0.85em', color: '#666', fontWeight: 500 }}>
                                Min QUAL:
                            </label>
                            <input 
                                type="number" 
                                step="10" 
                                min="0" 
                                value={minQualThreshold} 
                                onChange={(e) => setMinQualThreshold(e.target.value)}
                                title="Minimum QUAL score (Phred-scaled confidence) to include a GATK variant. 0 = no filter. Common thresholds: 30 (moderate), 100 (high confidence)"
                                style={{ 
                                    padding: '6px', 
                                    borderRadius: '4px', 
                                    border: '1px solid #ccc', 
                                    width: '70px',
                                    textAlign: 'center'
                                }} 
                            />
                        </div>
                    )}
                </div>
            </div>

            <div style={{ display: 'flex', justifyContent: 'flex-end', gap: '10px', marginTop: '10px' }}>
                <button 
                    onClick={async () => {
                        try {
                            const response = await axios.get(`/api/v1/analysis/consensus/${jobId}`, {
                                params: { freq_threshold: parseFloat(consensusThreshold), min_qual: parseFloat(minQualThreshold) },
                                responseType: 'blob'
                            });
                            const url = window.URL.createObjectURL(new Blob([response.data]));
                            const link = document.createElement('a');
                            link.href = url;
                            link.setAttribute('download', `consensus_${jobId}.fasta`); 
                            document.body.appendChild(link);
                            link.click();
                            link.parentNode?.removeChild(link);
                        } catch (e) {
                            console.error("Failed to download consensus", e);
                            alert("Failed to generate consensus sequence. Make sure analysis is complete.");
                        }
                    }}
                    style={{
                        backgroundColor: '#2e7d32',
                        color: 'white',
                        border: 'none',
                        padding: '10px 20px',
                        borderRadius: '5px',
                        cursor: 'pointer',
                        fontSize: '1em'
                    }}
                >
                    Download Consensus FASTA
                </button>
                <button 
                    onClick={async () => {
                        try {
                            const response = await axios.get(`/api/v1/analysis/consensus/preview/${jobId}`, {
                                params: { freq_threshold: parseFloat(consensusThreshold), min_qual: parseFloat(minQualThreshold) }
                            });
                            setConsensusProposals(response.data);
                            setShowConsensusModal(true);
                        } catch (e) {
                            console.error("Failed to preview consensus", e);
                            alert("Failed to load consensus preview. Make sure analysis is complete.");
                        }
                    }}
                    style={{
                        backgroundColor: '#1976d2',
                        color: 'white',
                        border: 'none',
                        padding: '10px 20px',
                        borderRadius: '5px',
                        cursor: 'pointer',
                        fontSize: '1em'
                    }}
                >
                    Download Consensus GBK
                </button>
            </div>
            
            {results?.chromosome_copy_number && results.chromosome_copy_number.length > 0 && (
                <CopyNumberPanel
                    data={results.chromosome_copy_number}
                    aneuploidy={results.aneuploidy || []}
                    tool={results.cnv_tool}
                    displayName={getDisplayName}
                    selectedChrom={selectedChrom}
                    onSelectChrom={(seqId) => {
                        if (!chromosomes.includes(seqId)) return;
                        setSelectedChrom(seqId);
                        setSelectedVariantPos(null);
                        setViewCenterTop(null);
                        setViewCenterBottom(null);
                        setZoomWindowTop(5000);
                        setZoomWindowBottom(5000);
                        document.getElementById('coverage-section')?.scrollIntoView({ behavior: 'smooth', block: 'start' });
                    }}
                />
            )}

            {coverage.length > 0 && (
                <div id="coverage-section" style={{ marginBottom: '40px', padding: '20px', backgroundColor: '#fff', borderRadius: '8px', boxShadow: '0 2px 4px rgba(0,0,0,0.1)' }}>
                    <h3 style={{ color: '#2c3e50', marginBottom: '5px' }}>Genome Coverage</h3>
                    <p style={{ fontSize: '0.85em', color: '#888', margin: '0 0 10px 0' }}>
                        Click anywhere on the coverage plot to navigate the annotation browser below.
                    </p>
                    <div style={{ width: '100%', height: 400 }}>
                        <ResponsiveContainer width="100%" height="100%">
                            <ComposedChart 
                                data={chartCoverageData}
                                onClick={(e: any) => {
                                    if (e && e.activeLabel != null) {
                                        const clickedPos = typeof e.activeLabel === 'number' ? e.activeLabel : parseInt(e.activeLabel);
                                        if (!isNaN(clickedPos)) {
                                            setViewCenterBottom(clickedPos);
                                            setSelectedVariantPos(null);
                                        }
                                    }
                                }}
                                style={{ cursor: 'crosshair' }}
                            >
                                <CartesianGrid strokeDasharray="3 3" />
                                <XAxis 
                                    dataKey="pos" 
                                    type="number"
                                    domain={[0, maxPos]}
                                    label={{ value: 'Position (bp)', position: 'insideBottomRight', offset: -10 }} 
                                    tickFormatter={(value) => {
                                        if (value === 0) return '0';
                                        if (value < 1000000) return `${(value / 1000).toFixed(0)}k`;
                                        return `${(value / 1000000).toFixed(2)}M`;
                                    }}
                                />
                                <YAxis label={{ value: 'Depth', angle: -90, position: 'insideLeft' }} />
                                <Tooltip content={<CustomTooltip />} isAnimationActive={false} />
                                <Legend />
                                
                                {/* Viewport indicator showing annotation browser position */}
                                {viewCenterBottom !== null && (
                                    <ReferenceArea
                                        x1={zoomDomainBottom[0]}
                                        x2={zoomDomainBottom[1]}
                                        fill="#1976d2"
                                        fillOpacity={0.12}
                                        stroke="#1976d2"
                                        strokeOpacity={0.4}
                                        strokeWidth={1}
                                        ifOverflow="hidden"
                                    />
                                )}

                                {/* Coverage Area - data inherited from ComposedChart */}
                                <Area 
                                    type="step" 
                                    dataKey="depth" 
                                    stroke="#8884d8" 
                                    fill="#8884d8" 
                                    fillOpacity={0.3}
                                    name="Coverage" 
                                    isAnimationActive={false}
                                />
                            </ComposedChart>
                        </ResponsiveContainer>
                    </div>
                </div>
            )}
            
            {/* Separate Scatter Plot for Variants - independent of coverage */}
            {variantData.length > 0 && maxPos > 0 && (
                        <div style={{ marginTop: '30px' }}>
                            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '10px' }}>
                                <h4 style={{ color: '#2c3e50', margin: 0 }}>Detected Variants Distribution</h4>
                                {/* Zoom Controls for Top Plot */}
                                <div style={{ display: 'flex', gap: '5px' }}>
                                    <button
                                        onClick={() => setViewCenterTop(prev => (prev === null ? null : Math.max(zoomWindowTop / 2, prev - (zoomWindowTop * 0.2))))}
                                        title="Pan Left"
                                        style={{ padding: '5px 10px', cursor: 'pointer', fontWeight: 'bold' }}
                                    >
                                        {'<<'}
                                    </button>
                                    <button
                                        onClick={() => {
                                            const newWindow = Math.min(maxPos, zoomWindowTop * 1.5);
                                            setZoomWindowTop(newWindow);
                                            if (newWindow >= maxPos) {
                                                setViewCenterTop(maxPos / 2);
                                            } else {
                                                // Clamp center so the view stays within bounds
                                                setViewCenterTop(prev => {
                                                    if (prev === null) return null;
                                                    return Math.max(newWindow / 2, Math.min(maxPos - newWindow / 2, prev));
                                                });
                                            }
                                        }}
                                        title="Zoom Out"
                                        style={{ padding: '5px 10px', cursor: 'pointer', fontWeight: 'bold' }}
                                    >
                                        -
                                    </button>
                                    <span style={{ padding: '5px', fontSize: '0.9em', alignSelf: 'center', color: '#666', minWidth: '80px', textAlign: 'center' }}>
                                        {zoomWindowTop < 1000 ? `${zoomWindowTop} bp` : `${(zoomWindowTop/1000).toFixed(1)} kb`}
                                    </span>
                                    <button
                                        onClick={() => setZoomWindowTop(curr => Math.max(100, curr / 1.5))}
                                        title="Zoom In"
                                        style={{ padding: '5px 10px', cursor: 'pointer', fontWeight: 'bold' }}
                                    >
                                        +
                                    </button>
                                    <button
                                        onClick={() => setViewCenterTop(prev => (prev === null ? null : Math.min(maxPos - zoomWindowTop / 2, prev + (zoomWindowTop * 0.2))))}
                                        title="Pan Right"
                                        style={{ padding: '5px 10px', cursor: 'pointer', fontWeight: 'bold' }}
                                    >
                                        {'>>'}
                                    </button>
                                </div>
                            </div>

                            <p style={{ fontSize: '0.9em', color: '#666' }}>
                                Click on a variant to center the Genome Annotation (bottom plot) on it.
                            </p>
                            
                            {/* Legend for Variant Types */}
                            <div style={{ display: 'flex', gap: '20px', marginBottom: '10px', fontSize: '0.9em' }}>
                                <div style={{ display: 'flex', alignItems: 'center', gap: '5px' }}>
                                    <div style={{ width: '10px', height: '10px', borderRadius: '50%', backgroundColor: '#d32f2f' }}></div>
                                    <span>SNP</span>
                                </div>
                                <div style={{ display: 'flex', alignItems: 'center', gap: '5px' }}>
                                    <div style={{ width: '10px', height: '10px', borderRadius: '50%', backgroundColor: '#1976d2' }}></div>
                                    <span>Deletion (DEL)</span>
                                </div>
                                <div style={{ display: 'flex', alignItems: 'center', gap: '5px' }}>
                                    <div style={{ width: '10px', height: '10px', borderRadius: '50%', backgroundColor: '#388e3c' }}></div>
                                    <span>Insertion (INS)</span>
                                </div>
                                <div style={{ display: 'flex', alignItems: 'center', gap: '5px' }}>
                                    <div style={{ width: '10px', height: '10px', borderRadius: '50%', backgroundColor: '#fbc02d' }}></div>
                                    <span>Substitution (SUB)</span>
                                </div>
                            </div>

                            <div style={{ width: '100%', height: 200 }}>
                                <ResponsiveContainer width="100%" height="100%">
                                    <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 20 }}>
                                        <CartesianGrid />
                                        <XAxis 
                                            type="number" 
                                            dataKey="x" 
                                            name="Position" 
                                            domain={zoomDomainTop}
                                            tickFormatter={(value) => {
                                                if (value === 0) return '0';
                                                if (value < 1000000) return `${(value / 1000).toFixed(0)}k`;
                                                return `${(value / 1000000).toFixed(2)}M`;
                                            }}
                                            allowDataOverflow={true}
                                            label={{ value: 'Position (bp)', position: 'insideBottomRight', offset: -10 }}
                                        />
                                        <YAxis 
                                            type="number" 
                                            dataKey="y" 
                                            name="Frequency" 
                                            unit="%" 
                                            domain={[0, 100]}
                                            label={{ value: 'Freq (%)', angle: -90, position: 'insideLeft' }}
                                        />
                                        <Tooltip cursor={{ strokeDasharray: '3 3' }} content={<CustomTooltip />} isAnimationActive={false} />
                                        <Scatter 
                                            name="Variants" 
                                            data={variantData} 
                                            shape={<VariantShape />}
                                            onClick={(data: any) => {
                                                if (data && data.payload && data.payload.pos) {
                                                    setSelectedVariantPos(data.payload.pos);
                                                } else if (data && data.pos) {
                                                    // Fallback for direct data object
                                                    setSelectedVariantPos(data.pos);
                                                }
                                            }}
                                            style={{ cursor: 'pointer' }}
                                        >
                                            {variantData.map((entry, index) => (
                                                <Cell 
                                                    key={`cell-${index}`} 
                                                    fill={entry.fill} 
                                                    stroke={entry.pos === selectedVariantPos ? '#000' : 'none'}
                                                    strokeWidth={entry.pos === selectedVariantPos ? 2 : 0}
                                                />
                                            ))}
                                        </Scatter>
                                    </ScatterChart>
                                </ResponsiveContainer>
                            </div>
                        </div>
                    )}

            {/* Genome Features Track - SVG-based, independent of coverage */}
            {(featureData.length > 0 || variantData.length > 0) && maxPos > 0 && (
                <>
                    <GenomeBrowser
                        features={features.filter(f => !selectedChrom || f.seq_id === selectedChrom || !f.seq_id)}
                        variants={variantData}
                        viewStart={zoomDomainBottom[0]}
                        viewEnd={zoomDomainBottom[1]}
                        maxPos={maxPos}
                        onViewChange={(start, end) => {
                            const center = (start + end) / 2;
                            const window = end - start;
                            setViewCenterBottom(center);
                            setZoomWindowBottom(window);
                        }}
                        onFeatureClick={openEditModal}
                        onVariantClick={(v) => {
                            setSelectedVariantPos(v.pos);
                        }}
                        selectedVariantPos={selectedVariantPos}
                        jobId={jobId}
                        onEvidenceOpen={(url) => setEvidenceUrl(
                            results?.variant_caller === 'gatk'
                                ? `/api/v1/analysis/evidence/${jobId}?path=${encodeURIComponent(url)}`
                                : `/results_files/${jobId}/breseq/output/${url}`
                        )}
                    />

                    {/* Selected Variant Details & Actions */}
                    {selectedVariantPos !== null && (
                        <div style={{ 
                            marginTop: '20px', 
                            padding: '15px', 
                            borderRadius: '8px', 
                            border: '1px solid #e0e0e0',
                            backgroundColor: '#fafafa',
                            display: 'flex',
                            justifyContent: 'space-between',
                            alignItems: 'center'
                        }}>
                            <div>
                                <h5 style={{ margin: '0 0 5px 0', color: '#1565c0' }}>
                                    Selected Position: {selectedVariantPos.toLocaleString()} bp
                                </h5>
                                {(() => {
                                    // Find features overlapping this position
                                    const overlapping = features.filter(f => 
                                        (!selectedChrom || f.seq_id === selectedChrom || !f.seq_id) &&
                                        selectedVariantPos >= f.start && selectedVariantPos <= f.end
                                    );
                                    
                                    if (overlapping.length > 0) {
                                        return (
                                            <div>
                                                <span style={{ fontWeight: 'bold', color: '#2e7d32' }}>Overlapping Gene(s): </span>
                                                {overlapping.map((f, i) => (
                                                    <div key={i} style={{ display: 'inline-block', marginRight: '15px' }}>
                                                        <span style={{ fontWeight: 500 }}>{f.name || f.locus_tag}</span>
                                                        <span style={{ fontSize: '0.9em', color: '#666', marginLeft: '5px' }}>({f.product})</span>
                                                        <button 
                                                            onClick={() => openEditModal(f)}
                                                            style={{ 
                                                                marginLeft: '8px',
                                                                padding: '2px 8px',
                                                                fontSize: '0.8em',
                                                                border: '1px solid #1976d2',
                                                                borderRadius: '4px',
                                                                background: 'white',
                                                                color: '#1976d2',
                                                                cursor: 'pointer'
                                                            }}
                                                        >
                                                            Edit
                                                        </button>
                                                    </div>
                                                ))}
                                            </div>
                                        );
                                    } else {
                                        return (
                                            <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
                                                <span style={{ color: '#666', fontSize: '0.9em' }}>No annotation overlap.</span>
                                                <button 
                                                    onClick={() => {
                                                        setAddForm({
                                                            name: '', locus_tag: '', product: '',
                                                            start: (selectedVariantPos - 50).toString(),
                                                            end: (selectedVariantPos + 50).toString(),
                                                            seq_id: selectedChrom || chromosomes[0] || 'unknown',
                                                            strand: '1'
                                                        });
                                                        setIsAddModalOpen(true);
                                                    }}
                                                    style={{ 
                                                        padding: '6px 12px',
                                                        fontSize: '0.9em',
                                                        backgroundColor: '#1976d2',
                                                        color: 'white',
                                                        border: 'none',
                                                        borderRadius: '4px',
                                                        cursor: 'pointer',
                                                        display: 'flex',
                                                        alignItems: 'center',
                                                        gap: '5px'
                                                    }}
                                                >
                                                    <span>+</span> Add Annotation
                                                </button>
                                            </div>
                                        );
                                    }
                                })()}
                            </div>
                            
                            {/* Show variant stats if a variant exists exactly at this pos */}
                            {(() => {
                                const exactVar = variantData.find(v => v.pos === selectedVariantPos);
                                if (exactVar) {
                                    return (
                                        <div style={{ textAlign: 'right', fontSize: '0.9em', marginLeft: '20px' }}>
                                            <div style={{ fontWeight: 'bold', color: '#1565c0' }}>{exactVar.mutation}</div>
                                            <div>Freq: <span style={{fontWeight: 'bold'}}>{exactVar.freq}%</span></div>
                                            <div style={{ color: '#666' }}>{exactVar.type}</div>
                                        </div>
                                    );
                                }
                                return null;
                            })()}
                        </div>
                    )}
                </>
            )}
            
            {results?.read_stats && (
                <div style={{ marginBottom: '30px' }}>
                    <h3 style={{ color: '#e65100' }}>Mapping Statistics</h3>
                    <div style={{ 
                        padding: '10px', 
                        backgroundColor: '#fff3e0', 
                        borderRadius: '4px',
                        border: '1px solid #ffe0b2'
                    }}>
                        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))', gap: '20px' }}>
                            <div>
                                <div style={{ fontSize: '0.9em', color: '#666' }}>Total Reads</div>
                                <div style={{ fontSize: '1.2em', fontWeight: 'bold' }}>{results.read_stats.total_reads.toLocaleString()}</div>
                            </div>
                            <div>
                                <div style={{ fontSize: '0.9em', color: '#666' }}>Mapped Reads</div>
                                <div style={{ fontSize: '1.2em', fontWeight: 'bold', color: '#4caf50' }}>{results.read_stats.mapped_percentage}%</div>
                            </div>
                            <div>
                                <div style={{ fontSize: '0.9em', color: '#666' }}>Unmapped Reads</div>
                                <div style={{ fontSize: '1.2em', fontWeight: 'bold', color: '#f44336' }}>
                                    {results.read_stats.unmapped_percentage}% (~{results.read_stats.unmapped_reads.toLocaleString()})
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            )}

            {results?.unmapped_files && results.unmapped_files.length > 0 && (
                <div style={{ marginBottom: '30px' }}>
                    <h3 style={{ color: '#e65100' }}>Unmapped Reads</h3>
                    
                    <ul style={{ listStyleType: 'none', padding: 0 }}>
                        {results.unmapped_files.map((file, idx) => (
                            <li key={idx} style={{ marginBottom: '8px' }}>
                                <a 
                                    href={`/results_files/${jobId}/${results.variant_caller === 'gatk' ? 'gatk' : 'breseq/data'}/${file}`}
                                    target="_blank"
                                    rel="noopener noreferrer"
                                    style={{ color: '#1976d2', textDecoration: 'none', display: 'flex', alignItems: 'center', gap: '8px' }}
                                >
                                    <span style={{ fontSize: '1.2em' }}>📄</span>
                                    {file}
                                </a>
                            </li>
                        ))}
                    </ul>
                </div>
            )}
            
            {results?.blast_results && results.blast_results.length > 0 && (
                <div style={{ marginBottom: '30px' }}>
                    <h3 style={{ color: '#7b1fa2' }}>Unmapped Reads Classification (BLAST Sample)</h3>
                    <p style={{ fontSize: '0.9em', color: '#666' }}>
                        Top hits from a random sample of unmapped reads against the NCBI Nucleotide database.
                    </p>
                    <table style={{ width: '100%', borderCollapse: 'collapse', marginTop: '10px', fontSize: '0.9em' }}>
                        <thead>
                            <tr style={{ backgroundColor: '#f2f2f2' }}>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Organism / Title</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Read Count</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Avg Identity</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Best E-value</th>
                            </tr>
                        </thead>
                        <tbody>
                            {results.blast_results.map((hit, idx) => (
                                <tr key={idx}>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{hit.title}</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{hit.count}</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{hit.identity}%</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{hit.evalue}</td>
                                </tr>
                            ))}
                        </tbody>
                    </table>
                </div>
            )}

            {results?.mash_results && results.mash_results.length > 0 && (
                <div style={{ marginBottom: '30px' }}>
                    <h3 style={{ color: '#d32f2f' }}>Contaminant Identification (MASH Screen)</h3>
                    <p style={{ fontSize: '0.9em', color: '#666' }}>
                        Screening input reads against RefSeq genomes to identify contained strains.
                    </p>
                    <table style={{ width: '100%', borderCollapse: 'collapse', marginTop: '10px', fontSize: '0.9em' }}>
                        <thead>
                            <tr style={{ backgroundColor: '#f2f2f2' }}>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Identity</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Shared Hashes</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Median Multiplicity</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>P-value</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Organism</th>
                            </tr>
                        </thead>
                        <tbody>
                            {results.mash_results.map((hit, idx) => (
                                <tr key={idx}>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{(hit.identity * 100).toFixed(2)}%</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{hit.shared_hashes}</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{hit.median_multiplicity}</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{hit.p_value}</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{hit.comment}</td>
                                </tr>
                            ))}
                        </tbody>
                    </table>
                </div>
            )}

            {results?.mash_unmapped_results && results.mash_unmapped_results.length > 0 && (
                <div style={{ marginBottom: '30px' }}>
                    <h3 style={{ color: '#d32f2f' }}>Contaminants in Unmapped Reads (MASH Screen)</h3>
                    <p style={{ fontSize: '0.9em', color: '#666' }}>
                        Screening only reads that did NOT map to the reference genome. This helps identify distant contaminants like S. aureus.
                    </p>
                    <table style={{ width: '100%', borderCollapse: 'collapse', marginTop: '10px', fontSize: '0.9em' }}>
                        <thead>
                            <tr style={{ backgroundColor: '#f2f2f2' }}>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Identity</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Shared Hashes</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Median Multiplicity</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>P-value</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Organism</th>
                            </tr>
                        </thead>
                        <tbody>
                            {results.mash_unmapped_results.map((hit, idx) => (
                                <tr key={idx}>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{(hit.identity * 100).toFixed(2)}%</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{hit.shared_hashes}</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{hit.median_multiplicity}</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{hit.p_value}</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{hit.comment}</td>
                                </tr>
                            ))}
                        </tbody>
                    </table>
                </div>
            )}

            {results?.contaminants && results.contaminants.length > 0 && (
                <div style={{ marginBottom: '30px' }}>
                    <h3 style={{ color: '#d32f2f' }}>Detected Contaminants (Kraken2)</h3>
                    <p style={{ fontSize: '0.9em', color: '#666' }}>
                        Taxonomic classification of unmapped reads.
                    </p>
                    <table style={{ width: '100%', borderCollapse: 'collapse', marginTop: '10px' }}>
                        <thead>
                            <tr style={{ backgroundColor: '#f2f2f2' }}>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Organism</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Rank</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Abundance</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Reads</th>
                            </tr>
                        </thead>
                        <tbody>
                            {results.contaminants.map((c, idx) => (
                                <tr key={idx}>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{c.name}</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{c.rank}</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{c.percentage}%</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{c.reads}</td>
                                </tr>
                            ))}
                        </tbody>
                    </table>
                </div>
            )}

            {results?.mutations && results.mutations.length > 0 && (
                <div>
                    <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                        <h3 style={{ color: '#1976d2' }}>Detected Variants (Mixed Populations)</h3>
                        {results.hc_snps_file && (
                            <a 
                                href={`/results_files/${jobId}/${results.hc_snps_file}`}
                                target="_blank"
                                rel="noopener noreferrer"
                                style={{ 
                                    fontSize: '0.9em', 
                                    color: '#fff', 
                                    backgroundColor: '#1976d2', 
                                    padding: '6px 12px', 
                                    borderRadius: '4px', 
                                    textDecoration: 'none' 
                                }}
                            >
                                Download High-Confidence SNPs (CSV)
                            </a>
                        )}
                    </div>
                    <table style={{ width: '100%', borderCollapse: 'collapse', marginTop: '10px', fontSize: '0.9em' }}>
                        <thead>
                            <tr style={{ backgroundColor: '#f2f2f2' }}>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Evidence</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Position</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Mutation</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Freq</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Gene</th>
                                <th style={{ padding: '8px', border: '1px solid #ddd' }}>Description</th>
                            </tr>
                        </thead>
                        <tbody>
                            {[...results.mutations]
                                .sort((a, b) => {
                                    const freqA = parseFloat(a.freq.replace('%', '')) || 0;
                                    const freqB = parseFloat(b.freq.replace('%', '')) || 0;
                                    return freqB - freqA;
                                })
                                .map((m, idx) => (
                                <tr key={idx}>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>
                                        <div style={{ display: 'flex', flexDirection: 'column', gap: '4px' }}>
                                            <span style={{ 
                                                fontWeight: 'bold', 
                                                color: m.type === 'SNP' ? '#2e7d32' : '#c62828',
                                                padding: '2px 6px',
                                                borderRadius: '4px',
                                                backgroundColor: '#f0f0f0',
                                                display: 'inline-block',
                                                textAlign: 'center'
                                            }}>
                                                {m.type}
                                            </span>
                                            {m.evidence_links && m.evidence_links.map((link, i) => (
                                                <button 
                                                    key={i}
                                                    onClick={() => setEvidenceUrl(
                                                        results.variant_caller === 'gatk'
                                                            ? `/api/v1/analysis/evidence/${jobId}?path=${encodeURIComponent(link.url)}`
                                                            : `/results_files/${jobId}/breseq/output/${link.url}`
                                                    )}
                                                    style={{ 
                                                        fontSize: '0.8em', 
                                                        color: '#1976d2', 
                                                        textDecoration: 'underline', 
                                                        border: 'none', 
                                                        background: 'none', 
                                                        padding: 0, 
                                                        cursor: 'pointer',
                                                        textAlign: 'left'
                                                    }}
                                                >
                                                    View Evidence
                                                </button>
                                            ))}
                                        </div>
                                    </td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{m.position}</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{m.mutation}</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>{m.freq}</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd', fontStyle: 'italic' }}>{m.gene}</td>
                                    <td style={{ padding: '8px', border: '1px solid #ddd' }}>
                                        <div>{m.description}</div>
                                        <div style={{ fontSize: '0.85em', color: '#666', marginTop: '4px' }}>{m.annotation}</div>
                                    </td>
                                </tr>
                            ))}
                        </tbody>
                    </table>
                </div>
            )}
            
            {(!results?.contaminants?.length && !results?.mutations?.length) && (
                <p>No contamination or variants detected.</p>
            )}

            {/* Evidence Modal */}
            {evidenceUrl && (
                <div style={{
                    position: 'fixed',
                    top: 0,
                    left: 0,
                    width: '100%',
                    height: '100%',
                    backgroundColor: 'rgba(0, 0, 0, 0.5)',
                    display: 'flex',
                    flexDirection: 'column',
                    justifyContent: 'center',
                    alignItems: 'center',
                    zIndex: 2000
                }}>
                    <div style={{
                        width: '90%',
                        height: '90%',
                        backgroundColor: '#fff',
                        borderRadius: '8px',
                        boxShadow: '0 4px 6px rgba(0,0,0,0.1)',
                        display: 'flex',
                        flexDirection: 'column',
                        overflow: 'hidden'
                    }}>
                        <div style={{
                            padding: '10px 20px',
                            borderBottom: '1px solid #eee',
                            display: 'flex',
                            justifyContent: 'space-between',
                            alignItems: 'center',
                            backgroundColor: '#f8f9fa'
                        }}>
                            <h3 style={{ margin: 0, fontSize: '1.1em', color: '#333' }}>Mutation Detailed Evidence</h3>
                            <button 
                                onClick={() => setEvidenceUrl(null)}
                                style={{
                                    border: 'none',
                                    background: 'none',
                                    fontSize: '1.5em',
                                    cursor: 'pointer',
                                    color: '#666'
                                }}
                            >
                                &times;
                            </button>
                        </div>
                        <div style={{ flex: 1, position: 'relative' }}>
                            <iframe 
                                src={evidenceUrl}
                                title="Evidence Viewer"
                                style={{ width: '100%', height: '100%', border: 'none' }}
                                onLoad={(e) => {
                                    // Inject styles into the iframe to match application theme
                                    const iframe = e.currentTarget;
                                    try {
                                        const doc = iframe.contentDocument || iframe.contentWindow?.document;
                                        if (doc) {
                                            const style = doc.createElement('style');
                                            style.textContent = `
                                                body {
                                                    font-family: system-ui, Avenir, Helvetica, Arial, sans-serif !important;
                                                    color: #333 !important;
                                                    background-color: #fff !important;
                                                    margin: 20px !important;
                                                }
                                                a { color: #1976d2 !important; text-decoration: none !important; }
                                                a:hover { text-decoration: underline !important; }
                                                table { border-collapse: collapse !important; width: 100% !important; margin-bottom: 20px !important; }
                                                th, td { 
                                                    padding: 12px 15px !important; 
                                                    border: 1px solid #e0e0e0 !important; 
                                                    text-align: left !important; 
                                                }
                                                th { background-color: #f5f5f5 !important; font-weight: 600 !important; color: #2c3e50 !important; }
                                                tr:nth-child(even) { background-color: #fafafa !important; }
                                                h1, h2, h3 { color: #2c3e50 !important; font-weight: 500 !important; }
                                                .polymorphism_alignment { font-family: monospace !important; background: #f8f9fa !important; padding: 10px !important; border-radius: 4px !important; overflow-x: auto !important; }
                                            `;
                                            doc.head.appendChild(style);
                                        }
                                    } catch (err) {
                                        console.warn("Could not inject styles into iframe (likely CORS restricted if separate domain, but should work for same-origin proxy)", err);
                                    }
                                }}
                            />
                        </div>
                    </div>
                </div>
            )}

            {/* Edit Feature Modal */}
            {editingFeature && (
                <div style={{
                    position: 'fixed', top: 0, left: 0, right: 0, bottom: 0,
                    backgroundColor: 'rgba(0,0,0,0.5)', 
                    display: 'flex', justifyContent: 'center', alignItems: 'center', zIndex: 2000
                }}>
                    <div style={{ backgroundColor: 'white', padding: '20px', borderRadius: '8px', minWidth: '400px', boxShadow: '0 4px 6px rgba(0,0,0,0.1)' }}>
                        <h3 style={{ marginTop: 0 }}>Edit Annotation</h3>
                        <p style={{fontSize: '0.9em', color: '#666'}}>
                            Adjusting coordinates is not supported for existing genes. Delete and add a new one if needed.
                        </p>
                        <div style={{ marginBottom: '15px' }}>
                            <label style={{display: 'block', marginBottom: '5px', fontWeight: 'bold'}}>Gene Name</label>
                            <input 
                                style={{width: '100%', padding: '8px', borderRadius: '4px', border: '1px solid #ddd'}}
                                value={editForm.name} 
                                onChange={e => setEditForm({...editForm, name: e.target.value})} 
                            />
                        </div>
                        <div style={{ marginBottom: '15px' }}>
                            <label style={{display: 'block', marginBottom: '5px', fontWeight: 'bold'}}>Locus Tag</label>
                            <input 
                                style={{width: '100%', padding: '8px', borderRadius: '4px', border: '1px solid #ddd'}}
                                value={editForm.locus_tag} 
                                onChange={e => setEditForm({...editForm, locus_tag: e.target.value})} 
                            />
                        </div>
                        <div style={{ marginBottom: '15px' }}>
                            <label style={{display: 'block', marginBottom: '5px', fontWeight: 'bold'}}>Product</label>
                            <input 
                                style={{width: '100%', padding: '8px', borderRadius: '4px', border: '1px solid #ddd'}}
                                value={editForm.product} 
                                onChange={e => setEditForm({...editForm, product: e.target.value})} 
                            />
                        </div>
                        <div style={{ display: 'flex', justifyContent: 'flex-end', gap: '10px', marginTop: '20px' }}>
                            <button 
                                onClick={handleDeleteFeature}
                                style={{ padding: '8px 16px', backgroundColor: '#e53935', color: 'white', border: 'none', borderRadius: '4px', cursor: 'pointer', marginRight: 'auto' }}
                            >
                                Delete
                            </button>
                            <button 
                                onClick={closeEditModal}
                                style={{ padding: '8px 16px', backgroundColor: '#f5f5f5', border: '1px solid #ddd', borderRadius: '4px', cursor: 'pointer' }}
                            >
                                Cancel
                            </button>
                            <button 
                                onClick={handleSaveEdit}
                                style={{ padding: '8px 16px', backgroundColor: '#1976d2', color: 'white', border: 'none', borderRadius: '4px', cursor: 'pointer' }}
                            >
                                Save Changes
                            </button>
                        </div>
                    </div>
                </div>
            )}

            {/* Add Feature Modal */}
            {isAddModalOpen && (
                <div style={{
                    position: 'fixed', top: 0, left: 0, right: 0, bottom: 0,
                    backgroundColor: 'rgba(0,0,0,0.5)', 
                    display: 'flex', justifyContent: 'center', alignItems: 'center', zIndex: 2000
                }}>
                    <div style={{ backgroundColor: 'white', padding: '20px', borderRadius: '8px', minWidth: '400px', maxHeight: '90vh', overflowY: 'auto', boxShadow: '0 4px 6px rgba(0,0,0,0.1)' }}>
                        <h3 style={{ marginTop: 0 }}>Add New Annotation</h3>
                        
                        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '15px', marginBottom: '15px' }}>
                            <div>
                                <label style={{display: 'block', marginBottom: '5px', fontWeight: 'bold'}}>Start (bp)</label>
                                <input 
                                    type="number"
                                    style={{width: '100%', padding: '8px', borderRadius: '4px', border: '1px solid #ddd'}}
                                    value={addForm.start} 
                                    onChange={e => setAddForm({...addForm, start: e.target.value})} 
                                />
                            </div>
                            <div>
                                <label style={{display: 'block', marginBottom: '5px', fontWeight: 'bold'}}>End (bp)</label>
                                <input 
                                    type="number"
                                    style={{width: '100%', padding: '8px', borderRadius: '4px', border: '1px solid #ddd'}}
                                    value={addForm.end} 
                                    onChange={e => setAddForm({...addForm, end: e.target.value})} 
                                />
                            </div>
                        </div>

                        <div style={{ marginBottom: '15px' }}>
                            <label style={{display: 'block', marginBottom: '5px', fontWeight: 'bold'}}>Strand</label>
                            <select 
                                style={{width: '100%', padding: '8px', borderRadius: '4px', border: '1px solid #ddd'}}
                                value={addForm.strand}
                                onChange={e => setAddForm({...addForm, strand: e.target.value})}
                            >
                                <option value="1">Forward (+)</option>
                                <option value="-1">Reverse (-)</option>
                            </select>
                        </div>

                        <div style={{ marginBottom: '15px' }}>
                            <label style={{display: 'block', marginBottom: '5px', fontWeight: 'bold'}}>Chromosome / Sequence ID</label>
                            <select 
                                style={{width: '100%', padding: '8px', borderRadius: '4px', border: '1px solid #ddd'}}
                                value={addForm.seq_id}
                                onChange={e => setAddForm({...addForm, seq_id: e.target.value})}
                            >
                                {chromosomes.length > 0 ? (
                                    chromosomes.map(chrom => (
                                        <option key={chrom} value={chrom}>{chrom}</option>
                                    ))
                                ) : (
                                    <option value="unknown">Unknown</option>
                                )}
                            </select>
                        </div>

                        <div style={{ marginBottom: '15px' }}>
                            <label style={{display: 'block', marginBottom: '5px', fontWeight: 'bold'}}>Gene Name</label>
                            <input 
                                style={{width: '100%', padding: '8px', borderRadius: '4px', border: '1px solid #ddd'}}
                                value={addForm.name} 
                                onChange={e => setAddForm({...addForm, name: e.target.value})} 
                            />
                        </div>
                        <div style={{ marginBottom: '15px' }}>
                            <label style={{display: 'block', marginBottom: '5px', fontWeight: 'bold'}}>Locus Tag</label>
                            <input 
                                style={{width: '100%', padding: '8px', borderRadius: '4px', border: '1px solid #ddd'}}
                                value={addForm.locus_tag} 
                                onChange={e => setAddForm({...addForm, locus_tag: e.target.value})} 
                            />
                        </div>
                        <div style={{ marginBottom: '15px' }}>
                            <label style={{display: 'block', marginBottom: '5px', fontWeight: 'bold'}}>Product</label>
                            <input 
                                style={{width: '100%', padding: '8px', borderRadius: '4px', border: '1px solid #ddd'}}
                                value={addForm.product} 
                                onChange={e => setAddForm({...addForm, product: e.target.value})} 
                            />
                        </div>
                        
                        <div style={{ display: 'flex', justifyContent: 'flex-end', gap: '10px', marginTop: '20px' }}>
                            <button 
                                onClick={() => setIsAddModalOpen(false)}
                                style={{ padding: '8px 16px', backgroundColor: '#f5f5f5', border: '1px solid #ddd', borderRadius: '4px', cursor: 'pointer' }}
                            >
                                Cancel
                            </button>
                            <button 
                                onClick={handleAddFeature}
                                style={{ padding: '8px 16px', backgroundColor: '#4caf50', color: 'white', border: 'none', borderRadius: '4px', cursor: 'pointer' }}
                            >
                                Add Feature
                            </button>
                        </div>
                    </div>
                </div>
            )}
            
            {showConsensusModal && (
                <div style={{
                    position: 'fixed', top: 0, left: 0, right: 0, bottom: 0,
                    backgroundColor: 'rgba(0,0,0,0.5)', 
                    display: 'flex', justifyContent: 'center', alignItems: 'center', zIndex: 2000
                }}>
                    <div style={{ backgroundColor: 'white', padding: '20px', borderRadius: '8px', width: '80%', maxHeight: '90vh', overflowY: 'auto', boxShadow: '0 4px 6px rgba(0,0,0,0.1)' }}>
                        <h3 style={{ marginTop: 0, marginBottom: '10px' }}>Consensus Genome Builder</h3>
                        <p style={{ fontSize: '0.9em', color: '#666', marginBottom: '20px' }}>
                            Select mutations to apply to annotation names. High confidence mutations (100%) are selected by default.
                            Sequence mutations will be applied based on threshold ({consensusThreshold}) regardless of name updates.
                        </p>
                        
                        <div style={{ maxHeight: '60vh', overflowY: 'auto', border: '1px solid #eee', marginBottom: '20px' }}>
                            <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: '0.9em' }}>
                                <thead style={{ backgroundColor: '#f5f5f5', position: 'sticky', top: 0 }}>
                                    <tr>
                                        <th style={{ padding: '10px', textAlign: 'left' }}>Select</th>
                                        <th style={{ padding: '10px', textAlign: 'left' }}>Freq</th>
                                        <th style={{ padding: '10px', textAlign: 'left' }}>Chr:Pos</th>
                                        <th style={{ padding: '10px', textAlign: 'left' }}>Context</th>
                                        <th style={{ padding: '10px', textAlign: 'left' }}>Type</th>
                                        <th style={{ padding: '10px', textAlign: 'left' }}>Original Name</th>
                                        <th style={{ padding: '10px', textAlign: 'left' }}>Proposed Name</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    {consensusProposals.map((p, idx) => (
                                        <tr key={idx} style={{ borderBottom: '1px solid #eee' }}>
                                            <td style={{ padding: '10px' }}>
                                                <input 
                                                    type="checkbox" 
                                                    checked={p.selected} 
                                                    onChange={e => {
                                                        const updated = [...consensusProposals];
                                                        updated[idx].selected = e.target.checked;
                                                        setConsensusProposals(updated);
                                                    }}
                                                />
                                            </td>
                                            <td style={{ padding: '10px' }}>{(p.frequency * 100).toFixed(1)}%</td>
                                            <td style={{ padding: '10px', whiteSpace: 'nowrap' }}>
                                                <div style={{ fontSize: '0.9em', color: '#666' }}>{getDisplayName(p.seq_id)}</div>
                                                <div style={{ fontWeight: 'bold' }}>{p.position}</div>
                                            </td>
                                            <td style={{ padding: '10px' }}>
                                                {p.nearby_features && p.nearby_features.length > 0 ? (
                                                    <div style={{ border: '1px solid #eee', backgroundColor: '#fafafa', padding: '2px', width: '200px' }}>
                                                        <svg width="200" height="40" style={{ display: 'block' }}>
                                                            {(() => {
                                                                // Use fixed window size +/- 2000 from center (pos)
                                                                const windowSize = 4000;
                                                                const center = p.position;
                                                                const start = center - 2000;
                                                                
                                                                return (
                                                                    <>
                                                                        {/* Axis Line */}
                                                                        <line x1="0" y1="20" x2="200" y2="20" stroke="#ccc" strokeWidth="1" />
                                                                        
                                                                        {/* Features */}
                                                                        {p.nearby_features.map((f: any, i: number) => {
                                                                            const scaleFn = (bp: number) => ((bp - start) / windowSize) * 200;
                                                                            // Calculate coordinates
                                                                            const x1_raw = scaleFn(f.start);
                                                                            const x2_raw = scaleFn(f.end);
                                                                            
                                                                            // Clamp to view
                                                                            const x1 = Math.max(0, x1_raw);
                                                                            const x2 = Math.min(200, x2_raw);
                                                                            // Ensure minimum visibility
                                                                            const width = Math.max(2, x2 - x1);
                                                                            
                                                                            // Determine if this is our target gene
                                                                            // The mutation is at center (p.position). Check if center is inside f.
                                                                            const isTarget = (p.position >= f.start && p.position <= f.end);
                                                                            
                                                                            const color = isTarget ? '#ff9800' : '#bbdefb'; // Orange for target, Blue for context
                                                                            const y = 14; 
                                                                            const h = 12;
                                                                            
                                                                            // Draw simple arrow shape
                                                                            // If gene is truncated by window, draw rect edge. If full end, draw arrow.
                                                                            // For simplicity, just draw direction if width allows
                                                                            let shape;
                                                                            
                                                                            if (width < 6) { // was 10
                                                                                shape = <rect x={x1} y={y} width={width} height={h} fill={color} stroke="#666" strokeWidth="0.5" opacity="0.8" />;
                                                                            } else {
                                                                                // Points for arrow
                                                                                // Forward (strand 1): >
                                                                                // Reverse (strand -1): <
                                                                                if (f.strand === 1) {
                                                                                    // Arrow points right
                                                                                    // If right edge is clipped (x2=200), don't draw point? Keep it simple.
                                                                                    shape = <polygon points={`${x1},${y} ${x2-5},${y} ${x2},${y+h/2} ${x2-5},${y+h} ${x1},${y+h}`} fill={color} stroke="#666" strokeWidth="0.5" opacity="0.8" />;
                                                                                } else {
                                                                                     // Arrow points left
                                                                                     shape = <polygon points={`${x1+5},${y} ${x2},${y} ${x2},${y+h} ${x1+5},${y+h} ${x1},${y+h/2}`} fill={color} stroke="#666" strokeWidth="0.5" opacity="0.8" />;
                                                                                }
                                                                            }

                                                                            return (
                                                                                <g key={i}>
                                                                                    {shape}
                                                                                    {/* Label if space permits */}
                                                                                    {width > 20 && ( // was 30
                                                                                        <text 
                                                                                            x={x1 + width/2} 
                                                                                            y={y-2} // Above the arrow
                                                                                            textAnchor="middle" 
                                                                                            fontSize="9" 
                                                                                            fontWeight="500"
                                                                                            fill="#444" 
                                                                                            style={{pointerEvents: 'none', userSelect: 'none'}}
                                                                                        >
                                                                                            {f.name}
                                                                                        </text>
                                                                                    )}
                                                                                </g>
                                                                            );
                                                                        })}
                                                                        
                                                                        {/* Mutation Marker (Red Dot) */}
                                                                        {/* Always at center (200 / 2 = 100) since we define window around it */}
                                                                         <circle cx="100" cy="10" r="3" fill="red" stroke="white" strokeWidth="1" />
                                                                         <line x1="100" y1="10" x2="100" y2="35" stroke="red" strokeWidth="1" strokeDasharray="2,2" opacity="0.5" />
                                                                    </>
                                                                );
                                                            })()}
                                                        </svg>
                                                    </div>
                                                ) : <span style={{ color: '#ccc', fontSize: '0.8em' }}>No Context</span>}
                                            </td>
                                            <td style={{ padding: '10px' }}>{p.is_coding ? 'Coding' : 'New'}</td>
                                            <td style={{ padding: '10px', color: '#666' }}>{p.gene_name || p.locus_tag || '-'}</td>
                                            <td style={{ padding: '10px' }}>
                                                <input
                                                    type="text"
                                                    value={p.proposed_name}
                                                    onChange={e => {
                                                        const updated = [...consensusProposals];
                                                        updated[idx].proposed_name = e.target.value;
                                                        // Auto-select if user edits the name manually
                                                        if (!updated[idx].selected) updated[idx].selected = true;
                                                        setConsensusProposals(updated);
                                                    }}
                                                    style={{ 
                                                        width: '100%', 
                                                        padding: '4px 8px',
                                                        border: '1px solid #ccc',
                                                        borderRadius: '4px',
                                                        fontWeight: 'bold', 
                                                        color: '#1565c0',
                                                        fontFamily: 'monospace'
                                                    }}
                                                />
                                            </td>
                                        </tr>
                                    ))}
                                    {consensusProposals.length === 0 && (
                                        <tr>
                                            <td colSpan={6} style={{ padding: '20px', textAlign: 'center', color: '#666' }}>
                                                No mutations found meeting threshold criteria.
                                            </td>
                                        </tr>
                                    )}
                                </tbody>
                            </table>
                        </div>

                        <div style={{ display: 'flex', justifyContent: 'flex-end', gap: '10px' }}>
                            <button 
                                onClick={() => setShowConsensusModal(false)}
                                style={{ padding: '10px 20px', backgroundColor: '#f5f5f5', border: '1px solid #ddd', borderRadius: '4px', cursor: 'pointer' }}
                            >
                                Cancel
                            </button>
                            <button 
                                onClick={async () => {
                                    try {
                                        const response = await axios.post(`/api/v1/analysis/consensus/download/${jobId}`, 
                                            consensusProposals, // Send as JSON body
                                            { responseType: 'blob', params: { freq_threshold: parseFloat(consensusThreshold), min_qual: parseFloat(minQualThreshold) } }
                                        );
                                        const url = window.URL.createObjectURL(new Blob([response.data]));
                                        const link = document.createElement('a');
                                        link.href = url;
                                        link.setAttribute('download', `consensus_${jobId}.gbk`); 
                                        document.body.appendChild(link);
                                        link.click();
                                        link.parentNode?.removeChild(link);
                                        setShowConsensusModal(false);
                                    } catch (e) {
                                        console.error("Failed to download consensus", e);
                                        alert("Failed to download consensus file.");
                                    }
                                }}
                                style={{ padding: '10px 20px', backgroundColor: '#1976d2', color: 'white', border: 'none', borderRadius: '4px', cursor: 'pointer' }}
                            >
                                Export / Download GBK
                            </button>
                        </div>
                    </div>
                </div>
            )}
        </div>
    );
};

export default ResultsView;
