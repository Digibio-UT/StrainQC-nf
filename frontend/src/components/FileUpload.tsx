import React, { useState, useEffect } from 'react';
import axios from 'axios';
import ResultsView from './ResultsView';

const FileUpload: React.FC = () => {
    // Reference Genome State
    const [refFile, setRefFile] = useState<File | null>(null);
    const [refUploading, setRefUploading] = useState(false);
    const [refFileId, setRefFileId] = useState<string>('');
    const [refMessage, setRefMessage] = useState<string>('');
    const [availableRefs, setAvailableRefs] = useState<{filename: string}[]>([]);
    const [useExistingRef, setUseExistingRef] = useState(true);

    // Sequencing Data State
    const [files, setFiles] = useState<FileList | null>(null);
    const [uploading, setUploading] = useState(false);
    const [message, setMessage] = useState<string>('');
    const [fileIds, setFileIds] = useState<string[]>([]);

    // Analysis State
    const [analysisMessage, setAnalysisMessage] = useState<string>('');
    const [jobId, setJobId] = useState<string>('');
    const [manualJobId, setManualJobId] = useState<string>('');
    
    // Analysis Options
    const [polymorphismMode, setPolymorphismMode] = useState<boolean>(false);
    const [variantCaller, setVariantCaller] = useState<'breseq' | 'gatk'>('breseq');
    const [ploidy, setPloidy] = useState<number>(1);
    const [showAdvanced, setShowAdvanced] = useState<boolean>(false);

    // breseq advanced options
    const [breseqLimitFoldCoverage, setBreseqLimitFoldCoverage] = useState<string>('');
    const [breseqBaseQualityCutoff, setBreseqBaseQualityCutoff] = useState<string>('3');
    const [breseqMinMappingQuality, setBreseqMinMappingQuality] = useState<string>('0');
    const [breseqPolyFreqCutoff, setBreseqPolyFreqCutoff] = useState<string>('0.05');
    const [breseqNoJunctionPrediction, setBreseqNoJunctionPrediction] = useState<boolean>(false);
    const [breseqConsensusCovEachStrand, setBreseqConsensusCovEachStrand] = useState<string>('0');
    const [breseqTargetedSequencing, setBreseqTargetedSequencing] = useState<boolean>(false);

    // GATK advanced options (shared between HC and Mutect2)
    const [gatkMinBaseQuality, setGatkMinBaseQuality] = useState<string>('10');
    const [gatkMinMappingQuality, setGatkMinMappingQuality] = useState<string>('20');
    const [gatkMaxReadsPerAlignStart, setGatkMaxReadsPerAlignStart] = useState<string>('50');
    const [gatkDontUseSoftClipped, setGatkDontUseSoftClipped] = useState<boolean>(false);
    // HC-specific
    const [gatkCallConf, setGatkCallConf] = useState<string>('30.0');
    // Mutect2-specific
    const [gatkCallableDepth, setGatkCallableDepth] = useState<string>('10');

    useEffect(() => {
        const fetchRefs = async () => {
            try {
                const res = await axios.get('/api/v1/upload/files?type=reference');
                setAvailableRefs(res.data);
            } catch (e) {
                console.error("Failed to fetch reference files", e);
            }
        };
        fetchRefs();
    }, []);

    const handleRefFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        if (e.target.files) setRefFile(e.target.files[0]);
    };

    const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        if (e.target.files) setFiles(e.target.files);
    };

    const handleLoadJob = () => {
        if (manualJobId) {
            setJobId(manualJobId);
            setAnalysisMessage(`Loaded job: ${manualJobId}`);
        }
    };

    const uploadFile = async (
        fileToUpload: File, 
        setLoading: (b: boolean) => void, 
        setMsg: (s: string) => void, 
        setId: (s: string) => void
    ) => {
        if (!fileToUpload) return;
        const formData = new FormData();
        formData.append('file', fileToUpload);

        setLoading(true);
        setMsg('');

        try {
            const response = await axios.post('/api/v1/upload/', formData, {
                headers: { 'Content-Type': 'multipart/form-data' },
            });
            setMsg(`Upload successful!`);
            setId(response.data.file_id);
        } catch (error) {
            setMsg('Upload failed. Please try again.');
            console.error(error);
        } finally {
            setLoading(false);
        }
    };

    const uploadFiles = async () => {
        if (!files || files.length === 0) return;
        setUploading(true);
        setMessage('');
        const newFileIds: string[] = [];
        
        try {
            for (let i = 0; i < files.length; i++) {
                const fileToUpload = files[i];
                const formData = new FormData();
                formData.append('file', fileToUpload);
                
                const response = await axios.post('/api/v1/upload/', formData, {
                    headers: { 'Content-Type': 'multipart/form-data' },
                });
                newFileIds.push(response.data.file_id);
            }
            setFileIds(newFileIds);
            setMessage(`Successfully uploaded ${newFileIds.length} files!`);
        } catch (error) {
            setMessage('Upload failed. Please try again.');
            console.error(error);
        } finally {
            setUploading(false);
        }
    };

    const handleStartAnalysis = async () => {
        if (fileIds.length === 0 || !refFileId) {
            setAnalysisMessage("Please upload both files first.");
            return;
        }
        try {
            const response = await axios.post('/api/v1/analysis/start', {
                file_ids: fileIds,
                reference_genome_id: refFileId,
                options: {
                    variant_caller: variantCaller,
                    polymorphism_mode: polymorphismMode,
                    ...(variantCaller === 'gatk' ? { ploidy } : {}),
                    // Advanced options
                    ...(variantCaller === 'breseq' ? {
                        breseq_limit_fold_coverage: breseqLimitFoldCoverage ? parseFloat(breseqLimitFoldCoverage) : 0,
                        breseq_base_quality_cutoff: parseInt(breseqBaseQualityCutoff) || 3,
                        breseq_min_mapping_quality: parseInt(breseqMinMappingQuality) || 0,
                        breseq_polymorphism_frequency_cutoff: parseFloat(breseqPolyFreqCutoff) || 0.05,
                        breseq_no_junction_prediction: breseqNoJunctionPrediction,
                        breseq_consensus_min_cov_each_strand: parseInt(breseqConsensusCovEachStrand) || 0,
                        breseq_targeted_sequencing: breseqTargetedSequencing,
                    } : {
                        gatk_min_base_quality: parseInt(gatkMinBaseQuality) || 10,
                        gatk_min_mapping_quality: parseInt(gatkMinMappingQuality) || 20,
                        gatk_max_reads_per_align_start: parseInt(gatkMaxReadsPerAlignStart) || 50,
                        gatk_dont_use_soft_clipped: gatkDontUseSoftClipped,
                        gatk_call_conf: parseFloat(gatkCallConf) || 30.0,
                        gatk_callable_depth: parseInt(gatkCallableDepth) || 10,
                    }),
                }
            });
            setAnalysisMessage(`Analysis started! Job ID: ${response.data.job_id}`);
            setJobId(response.data.job_id);
        } catch (error) {
            setAnalysisMessage('Failed to start analysis.');
            console.error(error);
        }
    }

    return (
        <div style={{ padding: '20px', border: '1px solid #ccc', borderRadius: '8px', width: '100%', maxWidth: '1200px', margin: '20px auto', textAlign: 'left' }}>
            
            {!jobId && (
                <>
                    {/* Load Existing Job */}
                    <div style={{ marginBottom: '20px', padding: '10px', backgroundColor: '#e3f2fd', borderRadius: '5px', display: 'flex', alignItems: 'center', gap: '10px' }}>
                        <h3 style={{ margin: 0, fontSize: '1em', whiteSpace: 'nowrap' }}>Load Existing Job:</h3>
                        <input 
                            type="text" 
                            placeholder="Enter Job ID (e.g., job_...)" 
                            value={manualJobId}
                            onChange={(e) => setManualJobId(e.target.value)}
                            style={{ flex: 1, padding: '5px' }}
                        />
                        <button onClick={handleLoadJob} disabled={!manualJobId}>Load</button>
                    </div>
                    
                    {/* Step 1: Reference Genome */}
                    <div style={{ marginBottom: '20px', padding: '10px', backgroundColor: '#f9f9f9', borderRadius: '5px' }}>
                        <h3>1. Reference Genome</h3>
                        <div style={{ marginBottom: '10px' }}>
                            <label style={{ marginRight: '15px' }}>
                                <input 
                                    type="radio" 
                                    checked={useExistingRef} 
                                    onChange={() => { setUseExistingRef(true); setRefFileId(''); setRefMessage(''); }} 
                                /> Select Existing
                            </label>
                            <label>
                                <input 
                                    type="radio" 
                                    checked={!useExistingRef} 
                                    onChange={() => { setUseExistingRef(false); setRefFileId(''); setRefMessage(''); }} 
                                /> Upload New
                            </label>
                        </div>

                        {useExistingRef ? (
                            <div>
                                <select 
                                    onChange={(e) => {
                                        setRefFileId(e.target.value);
                                        setRefMessage(e.target.value ? "Selected existing file" : "");
                                    }}
                                    value={refFileId}
                                    style={{ padding: '5px', width: '100%', maxWidth: '400px' }}
                                >
                                    <option value="">-- Select a reference genome --</option>
                                    {availableRefs.map(f => (
                                        <option key={f.filename} value={f.filename}>
                                            {f.filename}
                                        </option>
                                    ))}
                                </select>
                            </div>
                        ) : (
                            <div>
                                <p style={{ fontSize: '0.9em', color: '#666' }}>Accepted formats: .fasta, .fna, .gbk</p>
                                <input type="file" onChange={handleRefFileChange} accept=".fasta,.fna,.fa,.gbk,.gb" />
                                <button 
                                    onClick={() => refFile && uploadFile(refFile, setRefUploading, setRefMessage, setRefFileId)} 
                                    disabled={!refFile || refUploading}
                                    style={{ marginLeft: '10px' }}
                                >
                                    {refUploading ? 'Uploading...' : 'Upload Reference'}
                                </button>
                            </div>
                        )}
                        {refMessage && <span style={{ display: 'block', marginTop: '5px', color: refFileId ? 'green' : 'red' }}>{refMessage}</span>}
                    </div>

                    {/* Step 2: Sequencing Data */}
                    <div style={{ marginBottom: '20px', padding: '10px', backgroundColor: '#f9f9f9', borderRadius: '5px' }}>
                        <h3>2. Upload Sequencing Data (FASTQ)</h3>
                        <p style={{ fontSize: '0.9em', color: '#666' }}>Accepted formats: .fastq, .fastq.gz</p>
                        <input type="file" onChange={handleFileChange} accept=".fastq,.fastq.gz,.fq,.fq.gz,.gz" multiple />
                        <button 
                            onClick={uploadFiles} 
                            disabled={!files || uploading}
                            style={{ marginLeft: '10px' }}
                        >
                            {uploading ? 'Uploading...' : 'Upload Reads'}
                        </button>
                        {message && <span style={{ marginLeft: '10px', color: fileIds.length > 0 ? 'green' : 'red' }}>{message}</span>}
                    </div>

                    {/* Step 3: Analysis Pipeline Options — Tabbed */}
                    <div style={{ marginBottom: '20px', padding: '0', backgroundColor: '#f9f9f9', borderRadius: '5px', border: '1px solid #ddd' }}>
                        <h3 style={{ margin: '10px 10px 0 10px' }}>3. Pipeline Options</h3>
                        
                        {/* Tab Bar */}
                        <div style={{ display: 'flex', borderBottom: '2px solid #ddd', margin: '10px 10px 0 10px' }}>
                            <button
                                onClick={() => setVariantCaller('breseq')}
                                style={{
                                    padding: '10px 24px',
                                    cursor: 'pointer',
                                    border: 'none',
                                    borderBottom: variantCaller === 'breseq' ? '3px solid #1565c0' : '3px solid transparent',
                                    backgroundColor: variantCaller === 'breseq' ? '#e3f2fd' : 'transparent',
                                    color: variantCaller === 'breseq' ? '#1565c0' : '#666',
                                    fontWeight: variantCaller === 'breseq' ? 700 : 400,
                                    fontSize: '0.95em',
                                    borderRadius: '6px 6px 0 0',
                                    marginBottom: '-2px',
                                    transition: 'all 0.15s ease',
                                }}
                            >
                                breseq
                                <span style={{ fontSize: '0.8em', fontWeight: 400, marginLeft: '6px', color: '#888' }}>Microbial Evolution</span>
                            </button>
                            <button
                                onClick={() => setVariantCaller('gatk')}
                                style={{
                                    padding: '10px 24px',
                                    cursor: 'pointer',
                                    border: 'none',
                                    borderBottom: variantCaller === 'gatk' ? '3px solid #2e7d32' : '3px solid transparent',
                                    backgroundColor: variantCaller === 'gatk' ? '#e8f5e9' : 'transparent',
                                    color: variantCaller === 'gatk' ? '#2e7d32' : '#666',
                                    fontWeight: variantCaller === 'gatk' ? 700 : 400,
                                    fontSize: '0.95em',
                                    borderRadius: '6px 6px 0 0',
                                    marginBottom: '-2px',
                                    transition: 'all 0.15s ease',
                                }}
                            >
                                GATK
                                <span style={{ fontSize: '0.8em', fontWeight: 400, marginLeft: '6px', color: '#888' }}>BWA-MEM + HaplotypeCaller</span>
                            </button>
                        </div>

                        {/* Tab Content */}
                        <div style={{ padding: '15px' }}>
                            {variantCaller === 'breseq' ? (
                                /* ── breseq Tab ── */
                                <div>
                                    <p style={{ fontSize: '0.85em', color: '#555', margin: '0 0 12px 0' }}>
                                        Identifies mutations, structural variants, and new junctions from short-read resequencing data.
                                        Designed for haploid microbial genomes and evolved clones.
                                    </p>
                                    <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '8px' }}>
                                        <label style={{ cursor: 'pointer', display: 'flex', alignItems: 'center' }}>
                                            <input 
                                                type="checkbox" 
                                                checked={polymorphismMode} 
                                                onChange={(e) => setPolymorphismMode(e.target.checked)}
                                                style={{ marginRight: '8px' }}
                                            /> 
                                            Polymorphism Mode (Mixed Populations)
                                        </label>
                                        <span style={{ fontSize: '0.8em', color: '#888' }}>
                                            Uses the <code>-p</code> flag. Detects sub-clonal mutations. Requires more memory.
                                        </span>
                                    </div>
                                    {/* Advanced Options Toggle */}
                                    <div style={{ marginTop: '12px', borderTop: '1px solid #e0e0e0', paddingTop: '8px' }}>
                                        <button 
                                            onClick={() => setShowAdvanced(!showAdvanced)}
                                            style={{ background: 'none', border: 'none', cursor: 'pointer', color: '#1565c0', fontSize: '0.85em', padding: 0 }}
                                        >
                                            {showAdvanced ? '▾' : '▸'} Advanced Options
                                        </button>
                                        {showAdvanced && (
                                            <div style={{ marginTop: '10px', display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '10px', fontSize: '0.85em' }}>
                                                <label style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
                                                    <span>Limit Fold Coverage <span style={{ color: '#999' }}>(0 = unlimited)</span></span>
                                                    <input type="number" min="0" step="1" value={breseqLimitFoldCoverage} onChange={e => setBreseqLimitFoldCoverage(e.target.value)} placeholder="0" style={{ padding: '4px 6px', width: '120px' }} />
                                                </label>
                                                <label style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
                                                    <span>Base Quality Cutoff <span style={{ color: '#999' }}>(default: 3)</span></span>
                                                    <input type="number" min="0" max="40" value={breseqBaseQualityCutoff} onChange={e => setBreseqBaseQualityCutoff(e.target.value)} style={{ padding: '4px 6px', width: '120px' }} />
                                                </label>
                                                <label style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
                                                    <span>Min Mapping Quality <span style={{ color: '#999' }}>(default: 0)</span></span>
                                                    <input type="number" min="0" max="60" value={breseqMinMappingQuality} onChange={e => setBreseqMinMappingQuality(e.target.value)} style={{ padding: '4px 6px', width: '120px' }} />
                                                </label>
                                                <label style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
                                                    <span>Min Coverage Each Strand <span style={{ color: '#999' }}>(default: 0)</span></span>
                                                    <input type="number" min="0" value={breseqConsensusCovEachStrand} onChange={e => setBreseqConsensusCovEachStrand(e.target.value)} style={{ padding: '4px 6px', width: '120px' }} />
                                                </label>
                                                {polymorphismMode && (
                                                    <label style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
                                                        <span>Polymorphism Freq Cutoff <span style={{ color: '#999' }}>(default: 0.05)</span></span>
                                                        <input type="number" min="0" max="1" step="0.01" value={breseqPolyFreqCutoff} onChange={e => setBreseqPolyFreqCutoff(e.target.value)} style={{ padding: '4px 6px', width: '120px' }} />
                                                    </label>
                                                )}
                                                <label style={{ display: 'flex', alignItems: 'center', gap: '6px', cursor: 'pointer' }}>
                                                    <input type="checkbox" checked={breseqNoJunctionPrediction} onChange={e => setBreseqNoJunctionPrediction(e.target.checked)} />
                                                    <span>Skip Junction Prediction <span style={{ color: '#999' }}>(faster)</span></span>
                                                </label>
                                                <label style={{ display: 'flex', alignItems: 'center', gap: '6px', cursor: 'pointer' }}>
                                                    <input type="checkbox" checked={breseqTargetedSequencing} onChange={e => setBreseqTargetedSequencing(e.target.checked)} />
                                                    <span>Targeted Sequencing <span style={{ color: '#999' }}>(amplicon mode)</span></span>
                                                </label>
                                            </div>
                                        )}
                                    </div>
                                </div>
                            ) : (
                                /* ── GATK Tab ── */
                                <div>
                                    <p style={{ fontSize: '0.85em', color: '#555', margin: '0 0 12px 0' }}>
                                        Aligns reads with BWA-MEM and calls variants with GATK. 
                                        Suitable for both haploid and diploid organisms.
                                    </p>
                                    <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '12px' }}>
                                        <label style={{ cursor: 'pointer', display: 'flex', alignItems: 'center' }}>
                                            <input 
                                                type="checkbox" 
                                                checked={polymorphismMode} 
                                                onChange={(e) => setPolymorphismMode(e.target.checked)}
                                                style={{ marginRight: '8px' }}
                                            /> 
                                            Polymorphism Mode (Mixed Populations)
                                        </label>
                                        <span style={{ fontSize: '0.8em', color: '#888' }}>
                                            Uses Mutect2 instead of HaplotypeCaller to detect low-frequency variants.
                                        </span>
                                    </div>
                                    {!polymorphismMode && (
                                        <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
                                            <label style={{ cursor: 'pointer', display: 'flex', alignItems: 'center' }}>
                                                <strong style={{ marginRight: '8px' }}>Ploidy:</strong>
                                                <select 
                                                    value={ploidy} 
                                                    onChange={(e) => setPloidy(Number(e.target.value))}
                                                    style={{ padding: '4px 8px' }}
                                                >
                                                    <option value={1}>1 (Haploid)</option>
                                                    <option value={2}>2 (Diploid)</option>
                                                </select>
                                            </label>
                                        </div>
                                    )}
                                    {/* Advanced Options Toggle */}
                                    <div style={{ marginTop: '12px', borderTop: '1px solid #e0e0e0', paddingTop: '8px' }}>
                                        <button 
                                            onClick={() => setShowAdvanced(!showAdvanced)}
                                            style={{ background: 'none', border: 'none', cursor: 'pointer', color: '#2e7d32', fontSize: '0.85em', padding: 0 }}
                                        >
                                            {showAdvanced ? '▾' : '▸'} Advanced Options
                                        </button>
                                        {showAdvanced && (
                                            <div style={{ marginTop: '10px', display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '10px', fontSize: '0.85em' }}>
                                                <label style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
                                                    <span>Min Base Quality <span style={{ color: '#999' }}>(default: 10)</span></span>
                                                    <input type="number" min="0" max="40" value={gatkMinBaseQuality} onChange={e => setGatkMinBaseQuality(e.target.value)} style={{ padding: '4px 6px', width: '120px' }} />
                                                </label>
                                                <label style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
                                                    <span>Min Mapping Quality <span style={{ color: '#999' }}>(default: 20)</span></span>
                                                    <input type="number" min="0" max="60" value={gatkMinMappingQuality} onChange={e => setGatkMinMappingQuality(e.target.value)} style={{ padding: '4px 6px', width: '120px' }} />
                                                </label>
                                                <label style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
                                                    <span>Max Reads per Align Start <span style={{ color: '#999' }}>(default: 50)</span></span>
                                                    <input type="number" min="0" value={gatkMaxReadsPerAlignStart} onChange={e => setGatkMaxReadsPerAlignStart(e.target.value)} style={{ padding: '4px 6px', width: '120px' }} />
                                                </label>
                                                {!polymorphismMode && (
                                                    <label style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
                                                        <span>Call Confidence (QUAL) <span style={{ color: '#999' }}>(default: 30.0)</span></span>
                                                        <input type="number" min="0" step="1" value={gatkCallConf} onChange={e => setGatkCallConf(e.target.value)} style={{ padding: '4px 6px', width: '120px' }} />
                                                    </label>
                                                )}
                                                {polymorphismMode && (
                                                    <label style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
                                                        <span>Callable Depth <span style={{ color: '#999' }}>(default: 10)</span></span>
                                                        <input type="number" min="0" value={gatkCallableDepth} onChange={e => setGatkCallableDepth(e.target.value)} style={{ padding: '4px 6px', width: '120px' }} />
                                                    </label>
                                                )}
                                                <label style={{ display: 'flex', alignItems: 'center', gap: '6px', cursor: 'pointer' }}>
                                                    <input type="checkbox" checked={gatkDontUseSoftClipped} onChange={e => setGatkDontUseSoftClipped(e.target.checked)} />
                                                    <span>Don't Use Soft-Clipped Bases</span>
                                                </label>
                                            </div>
                                        )}
                                    </div>
                                </div>
                            )}
                        </div>
                    </div>

                    {/* Step 4: Action */}
                    <div style={{ borderTop: '1px solid #eee', paddingTop: '20px', textAlign: 'center' }}>
                        <button 
                            onClick={handleStartAnalysis}
                            disabled={fileIds.length === 0 || !refFileId || !!jobId}
                            style={{ padding: '10px 20px', fontSize: '1.1em', cursor: (fileIds.length === 0 || !refFileId || !!jobId) ? 'not-allowed' : 'pointer' }}
                        >
                            {jobId ? 'Analysis Running...' : 'Start analysis'}
                        </button>
                    </div>
                </>
            )}

            {analysisMessage && <p style={{ marginTop: '10px', fontWeight: 'bold', textAlign: 'center' }}>{analysisMessage}</p>}

            {/* Step 4: Results */}
            {jobId && (
                <div style={{ marginTop: '20px' }}>
                    <div style={{ marginBottom: '10px',  }}>
                        <a 
                            href="#" 
                            onClick={(e) => { 
                                e.preventDefault(); 
                                setJobId(''); 
                                setManualJobId('');
                                setAnalysisMessage('');
                                setFiles(null);
                                setRefFile(null);
                                setUseExistingRef(true);
                            }}
                            style={{ color: '#007bff', textDecoration: 'none', fontSize: '0.9em' }}
                        >
                            &larr; Back to Submission
                        </a>
                    </div>
                    <ResultsView jobId={jobId} />
                </div>
            )}
        </div>
    );
};

export default FileUpload;
