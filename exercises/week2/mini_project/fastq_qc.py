#!/usr/bin/env python3
"""FASTQ Quality Control (QC) Tool.

Reads a FASTQ file (gzip supported) and generates a premium one-page HTML report
with interactive charts for sequence lengths, quality scores, and base composition.
"""

import argparse
import os
import sys
import gzip
import json
from collections import Counter, defaultdict

# Standard genetic code dictionary not needed here since we are doing FASTQ QC

def generate_sample_fastq(filepath, num_reads=1000):
    """Generates a dummy FASTQ file for demonstration and testing."""
    import random
    
    bases = ['A', 'C', 'G', 'T']
    # A few potential adapters / overrepresented sequences
    adapters = ["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"]
    
    print(f"Generating dummy FASTQ file: {filepath} with {num_reads} reads...")
    with open(filepath, 'w') as f:
        for i in range(num_reads):
            # 95% normal reads, 5% adapter contaminated
            is_contaminated = random.random() < 0.05
            
            if is_contaminated:
                seq = random.choice(adapters) + "".join(random.choices(bases, k=random.randint(20, 50)))
            else:
                # normal sequence with varying quality
                length = random.randint(75, 100)
                seq = "".join(random.choices(bases, k=length))
            
            # Create a mock quality score sequence.
            # Normal reads start high quality and decrease slightly towards the end
            qual_scores = []
            for pos in range(len(seq)):
                # Base quality 30-40, slightly declining
                base_q = max(10, int(random.normalvariate(37 - (pos * 0.08), 3)))
                # clamp to 40
                base_q = min(40, base_q)
                qual_scores.append(chr(base_q + 33))
            
            qual_seq = "".join(qual_scores)
            
            f.write(f"@read_{i} mock_instrument:flowcell:lane:tile:x:y\n")
            f.write(f"{seq}\n")
            f.write(f"+\n")
            f.write(f"{qual_seq}\n")
    print("Dummy FASTQ file generation complete.")


def parse_fastq(filepath):
    """Parses FASTQ (supports gzip) and returns generator of (header, seq, qual)."""
    # Check if gzip
    if filepath.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
        
    with opener(filepath, mode) as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline().strip()
            yield header.strip(), seq, qual


def run_qc(fastq_path):
    """Computes quality control statistics for the FASTQ file."""
    total_reads = 0
    total_bases = 0
    gc_count = 0
    n_count = 0
    
    # Store distributions
    length_counter = Counter()
    quality_bins = Counter()  # Mean Q per read: <10, 10-19, 20-29, 30+
    
    # To track per-position stats, we need arrays up to the maximum read length
    max_read_len = 250  # Dynamic expansion
    per_pos_qual_sum = [0] * max_read_len
    per_pos_qual_cnt = [0] * max_read_len
    per_pos_base_counts = [defaultdict(int) for _ in range(max_read_len)]
    
    seq_counter = Counter()
    
    actual_max_len = 0
    
    print(f"Analyzing FASTQ file: {fastq_path}...")
    
    for header, seq, qual in parse_fastq(fastq_path):
        total_reads += 1
        seq_len = len(seq)
        total_bases += seq_len
        length_counter[seq_len] += 1
        
        if seq_len > actual_max_len:
            actual_max_len = seq_len
            # Expand lists if needed
            if seq_len > len(per_pos_qual_sum):
                diff = seq_len - len(per_pos_qual_sum)
                per_pos_qual_sum.extend([0] * diff)
                per_pos_qual_cnt.extend([0] * diff)
                per_pos_base_counts.extend([defaultdict(int) for _ in range(diff)])
        
        # Track sequence for overrepresentation (keep top 10000 to avoid memory blowup)
        if total_reads < 100000:
            seq_counter[seq] += 1
            
        read_qual_sum = 0
        for pos, (base, q_char) in enumerate(zip(seq, qual)):
            q_score = ord(q_char) - 33
            read_qual_sum += q_score
            
            per_pos_qual_sum[pos] += q_score
            per_pos_qual_cnt[pos] += 1
            
            base_upper = base.upper()
            per_pos_base_counts[pos][base_upper] += 1
            
            if base_upper == 'G' or base_upper == 'C':
                gc_count += 1
            elif base_upper == 'N':
                n_count += 1
                
        # Mean Q for this read
        if seq_len > 0:
            mean_q = read_qual_sum / seq_len
            if mean_q < 10:
                quality_bins['Q < 10'] += 1
            elif mean_q < 20:
                quality_bins['Q 10-19'] += 1
            elif mean_q < 30:
                quality_bins['Q 20-29'] += 1
            else:
                quality_bins['Q >= 30'] += 1
        else:
            quality_bins['Q < 10'] += 1
            
        if total_reads % 50000 == 0:
            print(f"  Processed {total_reads} reads...")
            
    if total_reads == 0:
        print("Error: No reads found in the FASTQ file.", file=sys.stderr)
        sys.exit(1)
        
    # Trim per-position arrays to actual max length
    per_pos_qual_sum = per_pos_qual_sum[:actual_max_len]
    per_pos_qual_cnt = per_pos_qual_cnt[:actual_max_len]
    per_pos_base_counts = per_pos_base_counts[:actual_max_len]
    
    # Calculate average quality per position
    per_pos_mean_qual = []
    for q_sum, q_cnt in zip(per_pos_qual_sum, per_pos_qual_cnt):
        per_pos_mean_qual.append(round(q_sum / q_cnt, 2) if q_cnt > 0 else 0)
        
    # Calculate base percentages per position
    per_pos_bases_pct = {'A': [], 'C': [], 'G': [], 'T': [], 'N': []}
    for pos_counts in per_pos_base_counts:
        pos_total = sum(pos_counts.values())
        for b in ['A', 'C', 'G', 'T', 'N']:
            pct = round((pos_counts[b] / pos_total) * 100, 2) if pos_total > 0 else 0
            per_pos_bases_pct[b].append(pct)
            
    # Calculate overall stats
    gc_content = round((gc_count / total_bases) * 100, 2) if total_bases > 0 else 0
    n_content = round((n_count / total_bases) * 100, 2) if total_bases > 0 else 0
    
    # Sequence lengths
    sorted_lengths = sorted(length_counter.keys())
    min_len = sorted_lengths[0]
    max_len = sorted_lengths[-1]
    avg_len = round(total_bases / total_reads, 1)
    
    # Overrepresented sequences (seqs with > 0.2% of total reads)
    overrepresented = []
    threshold = max(2, int(total_reads * 0.002))
    for seq, count in seq_counter.most_common(10):
        if count >= threshold:
            pct = round((count / total_reads) * 100, 3)
            overrepresented.append({
                'sequence': seq,
                'count': count,
                'percentage': pct
            })
            
    qc_stats = {
        'filename': os.path.basename(fastq_path),
        'total_reads': total_reads,
        'total_bases': total_bases,
        'gc_content': gc_content,
        'n_content': n_content,
        'min_len': min_len,
        'max_len': max_len,
        'avg_len': avg_len,
        'quality_bins': dict(quality_bins),
        'per_pos_mean_qual': per_pos_mean_qual,
        'per_pos_bases_pct': per_pos_bases_pct,
        'overrepresented': overrepresented,
        'length_distribution': {str(k): v for k, v in sorted(length_counter.items())}
    }
    
    return qc_stats


def generate_html_report(stats, output_path):
    """Generates a gorgeous standalone HTML report from the QC stats."""
    
    # Embed data as JSON inside the HTML template so Javascript can render interactive charts
    stats_json = json.dumps(stats)
    
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>FASTQ Quality Control Report - {stats['filename']}</title>
    <!-- Modern Premium Typography -->
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;500;600;700&family=Space+Grotesk:wght@400;500;600;700&display=swap" rel="stylesheet">
    
    <!-- Chart.js for premium interactive plots -->
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>

    <style>
        :root {{
            --bg-primary: #0b0f19;
            --bg-secondary: #161c2d;
            --bg-tertiary: #1f283e;
            --text-primary: #f3f4f6;
            --text-secondary: #9ca3af;
            --accent-blue: #3b82f6;
            --accent-cyan: #06b6d4;
            --accent-green: #10b981;
            --accent-red: #ef4444;
            --accent-purple: #8b5cf6;
            --border-color: rgba(255, 255, 255, 0.08);
            --card-shadow: 0 10px 30px -10px rgba(0, 0, 0, 0.5);
            --glow-blue: 0 0 20px rgba(59, 130, 246, 0.15);
        }}

        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}

        body {{
            font-family: 'Outfit', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            background-color: var(--bg-primary);
            color: var(--text-primary);
            line-height: 1.6;
            padding: 2rem;
            min-height: 100vh;
        }}

        .container {{
            max-width: 1400px;
            margin: 0 auto;
        }}

        header {{
            display: flex;
            align-items: center;
            justify-content: space-between;
            margin-bottom: 2.5rem;
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 1.5rem;
        }}

        .logo-title h1 {{
            font-family: 'Space Grotesk', sans-serif;
            font-size: 2.2rem;
            font-weight: 700;
            background: linear-gradient(135deg, var(--accent-cyan), var(--accent-blue));
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            margin-bottom: 0.25rem;
        }}

        .logo-title p {{
            color: var(--text-secondary);
            font-size: 0.95rem;
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }}

        .badge {{
            background: rgba(16, 185, 129, 0.15);
            color: var(--accent-green);
            padding: 0.25rem 0.75rem;
            border-radius: 50px;
            font-weight: 600;
            font-size: 0.8rem;
            border: 1px solid rgba(16, 185, 129, 0.3);
            letter-spacing: 0.05em;
            text-transform: uppercase;
        }}

        .meta-info {{
            font-size: 0.85rem;
            color: var(--text-secondary);
            text-align: right;
        }}

        /* Key Metrics Grid */
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 1.5rem;
            margin-bottom: 2.5rem;
        }}

        .metric-card {{
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            border-radius: 16px;
            padding: 1.5rem;
            box-shadow: var(--card-shadow);
            transition: transform 0.3s ease, border-color 0.3s ease;
            position: relative;
            overflow: hidden;
        }}

        .metric-card:hover {{
            transform: translateY(-4px);
            border-color: rgba(59, 130, 246, 0.3);
            box-shadow: var(--glow-blue);
        }}

        .metric-card::before {{
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            width: 4px;
            height: 100%;
            background: var(--accent-blue);
        }}

        .metric-card.green::before {{ background: var(--accent-green); }}
        .metric-card.purple::before {{ background: var(--accent-purple); }}
        .metric-card.cyan::before {{ background: var(--accent-cyan); }}

        .metric-label {{
            font-size: 0.85rem;
            color: var(--text-secondary);
            text-transform: uppercase;
            letter-spacing: 0.05em;
            margin-bottom: 0.5rem;
            font-weight: 500;
        }}

        .metric-value {{
            font-size: 1.8rem;
            font-weight: 700;
            font-family: 'Space Grotesk', sans-serif;
            color: var(--text-primary);
        }}

        .metric-sub {{
            font-size: 0.8rem;
            color: var(--text-secondary);
            margin-top: 0.25rem;
        }}

        /* Charts Layout */
        .charts-grid {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 1.5rem;
            margin-bottom: 2.5rem;
        }}

        @media (max-width: 1024px) {{
            .charts-grid {{
                grid-template-columns: 1fr;
            }}
            body {{
                padding: 1rem;
            }}
        }}

        .chart-card {{
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            border-radius: 16px;
            padding: 1.75rem;
            box-shadow: var(--card-shadow);
            display: flex;
            flex-direction: column;
        }}

        .chart-card-full {{
            grid-column: 1 / -1;
        }}

        .chart-header {{
            margin-bottom: 1.25rem;
        }}

        .chart-title {{
            font-size: 1.2rem;
            font-weight: 600;
            color: var(--text-primary);
            font-family: 'Space Grotesk', sans-serif;
        }}

        .chart-desc {{
            font-size: 0.85rem;
            color: var(--text-secondary);
            margin-top: 0.15rem;
        }}

        .chart-container {{
            position: relative;
            flex-grow: 1;
            min-height: 320px;
            display: flex;
            align-items: center;
            justify-content: center;
        }}

        /* Table Card */
        .table-card {{
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            border-radius: 16px;
            padding: 1.75rem;
            box-shadow: var(--card-shadow);
            margin-bottom: 2.5rem;
        }}

        .table-title {{
            font-size: 1.2rem;
            font-weight: 600;
            color: var(--text-primary);
            font-family: 'Space Grotesk', sans-serif;
            margin-bottom: 1rem;
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
            text-align: left;
        }}

        th {{
            border-bottom: 2px solid var(--border-color);
            padding: 0.75rem 1rem;
            font-weight: 600;
            font-size: 0.85rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: var(--text-secondary);
        }}

        td {{
            padding: 0.75rem 1rem;
            border-bottom: 1px solid var(--border-color);
            font-size: 0.9rem;
        }}

        tr:hover td {{
            background: rgba(255, 255, 255, 0.02);
            color: #fff;
        }}

        .mono {{
            font-family: 'Courier New', Courier, monospace;
            word-break: break-all;
            background: rgba(0, 0, 0, 0.2);
            padding: 0.15rem 0.4rem;
            border-radius: 4px;
            color: var(--accent-cyan);
        }}

        footer {{
            text-align: center;
            color: var(--text-secondary);
            font-size: 0.8rem;
            margin-top: 3rem;
            padding-top: 1.5rem;
            border-top: 1px solid var(--border-color);
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <div class="logo-title">
                <h1>FASTQ Quality Control</h1>
                <p>Report for file: <span style="font-weight:600; color:var(--text-primary);">{stats['filename']}</span> <span class="badge">Success</span></p>
            </div>
            <div class="meta-info">
                <p>Generated by Antigravity Bio QC</p>
                <p>Analysis Date: 2026-06-15</p>
            </div>
        </header>

        <!-- Metrics Grid -->
        <div class="metrics-grid">
            <div class="metric-card">
                <div class="metric-label">Total Reads</div>
                <div class="metric-value">{stats['total_reads']:,}</div>
                <div class="metric-sub">Sequenced clusters</div>
            </div>
            <div class="metric-card green">
                <div class="metric-label">Total Yield</div>
                <div class="metric-value">{stats['total_bases']:,}</div>
                <div class="metric-sub">Bases sequenced</div>
            </div>
            <div class="metric-card purple">
                <div class="metric-label">GC Content</div>
                <div class="metric-value">{stats['gc_content']}%</div>
                <div class="metric-sub">Target: ~40-60%</div>
            </div>
            <div class="metric-card cyan">
                <div class="metric-label">Mean Length</div>
                <div class="metric-value">{stats['avg_len']} bp</div>
                <div class="metric-sub">Range: {stats['min_len']} - {stats['max_len']} bp</div>
            </div>
            <div class="metric-card">
                <div class="metric-label">N Content</div>
                <div class="metric-value" style="color: { 'var(--accent-red)' if stats['n_content'] > 1 else 'var(--text-primary)' }">{stats['n_content']}%</div>
                <div class="metric-sub">Target: &lt; 1%</div>
            </div>
        </div>

        <!-- Charts Grid -->
        <div class="charts-grid">
            
            <!-- Quality Score Per Position -->
            <div class="chart-card chart-card-full">
                <div class="chart-header">
                    <div class="chart-title">Quality Scores Across Positions</div>
                    <div class="chart-desc">Mean Phred quality score (Q) at each cycle/position. Scores &gt; 30 are excellent (1 in 1000 error rate).</div>
                </div>
                <div class="chart-container">
                    <canvas id="qualityChart"></canvas>
                </div>
            </div>

            <!-- Base Composition Per Position -->
            <div class="chart-card">
                <div class="chart-header">
                    <div class="chart-title">Base Composition Per Cycle</div>
                    <div class="chart-desc">Percentage of A, C, G, T, and N bases at each cycle. A/T and G/C lines should stay parallel.</div>
                </div>
                <div class="chart-container">
                    <canvas id="basesChart"></canvas>
                </div>
            </div>

            <!-- Quality score distribution (read-level) -->
            <div class="chart-card">
                <div class="chart-header">
                    <div class="chart-title">Per-Read Mean Quality Distribution</div>
                    <div class="chart-desc">Percentage of reads falling into average Phred score bins. Ideally, most reads are Q30+.</div>
                </div>
                <div class="chart-container">
                    <canvas id="qualDistributionChart"></canvas>
                </div>
            </div>

            <!-- Read Length Distribution -->
            <div class="chart-card chart-card-full">
                <div class="chart-header">
                    <div class="chart-title">Sequence Length Distribution</div>
                    <div class="chart-desc">Frequency of read lengths. Peak indicates main fragment size range.</div>
                </div>
                <div class="chart-container">
                    <canvas id="lengthChart"></canvas>
                </div>
            </div>
            
        </div>

        <!-- Overrepresented Sequences -->
        <div class="table-card">
            <div class="table-title">Overrepresented Sequences (Top 10)</div>
            {f'''
            <table>
                <thead>
                    <tr>
                        <th style="width: 60%">Sequence</th>
                        <th>Count</th>
                        <th>Percentage</th>
                        <th>Status</th>
                    </tr>
                </thead>
                <tbody>
                    {"".join(f"""
                    <tr>
                        <td><span class="mono">{item['sequence']}</span></td>
                        <td>{item['count']:,}</td>
                        <td>{item['percentage']}%</td>
                        <td>
                            <span class="badge" style="background:rgba(239,68,68,0.15); color:var(--accent-red); border-color:rgba(239,68,68,0.3)">Overrepresented</span>
                        </td>
                    </tr>
                    """ for item in stats['overrepresented']) if stats['overrepresented'] else '<tr><td colspan="4" style="text-align:center; color:var(--text-secondary);">No overrepresented sequences found (all sequences are < 0.2% of total reads)</td></tr>'}
                </tbody>
            </table>
            '''}
        </div>

        <footer>
            <p>Antigravity Bioinformatics Suite • Version 2.0 • Brno Prep Log</p>
        </footer>
    </div>

    <!-- Data Injection for Charts -->
    <script>
        const qcData = {stats_json};
        
        // Chart.js Default Dark styling
        Chart.defaults.color = '#9ca3af';
        Chart.defaults.borderColor = 'rgba(255, 255, 255, 0.08)';
        Chart.defaults.font.family = "'Outfit', sans-serif";

        // 1. Quality Chart
        const cycles = Array.from({{length: qcData.per_pos_mean_qual.length}}, (_, i) => i + 1);
        new Chart(document.getElementById('qualityChart'), {{
            type: 'line',
            data: {{
                labels: cycles,
                datasets: [{{
                    label: 'Mean Quality Score (Q)',
                    data: qcData.per_pos_mean_qual,
                    borderColor: '#3b82f6',
                    backgroundColor: 'rgba(59, 130, 246, 0.1)',
                    fill: true,
                    tension: 0.2,
                    borderWidth: 2,
                    pointRadius: cycles.length > 100 ? 0 : 2
                }}]
            }},
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                scales: {{
                    y: {{
                        min: 0,
                        max: 42,
                        title: {{
                            display: true,
                            text: 'Phred Quality Score (Q)'
                        }}
                    }},
                    x: {{
                        title: {{
                            display: true,
                            text: 'Cycle / Position (bp)'
                        }}
                    }}
                }},
                plugins: {{
                    legend: {{ display: false }},
                    // Color background zones
                    annotation: {{
                        // Can be added later if needed, simple text annotations or colored stripes
                    }}
                }}
            }}
        }});

        // 2. Base Composition Chart
        new Chart(document.getElementById('basesChart'), {{
            type: 'line',
            data: {{
                labels: cycles,
                datasets: [
                    {{ label: 'A', data: qcData.per_pos_bases_pct.A, borderColor: '#10b981', borderWidth: 2, tension: 0.1, pointRadius: 0 }},
                    {{ label: 'T', data: qcData.per_pos_bases_pct.T, borderColor: '#ef4444', borderWidth: 2, tension: 0.1, pointRadius: 0 }},
                    {{ label: 'C', data: qcData.per_pos_bases_pct.C, borderColor: '#8b5cf6', borderWidth: 2, tension: 0.1, pointRadius: 0 }},
                    {{ label: 'G', data: qcData.per_pos_bases_pct.G, borderColor: '#3b82f6', borderWidth: 2, tension: 0.1, pointRadius: 0 }},
                    {{ label: 'N', data: qcData.per_pos_bases_pct.N, borderColor: '#e5e7eb', borderWidth: 2, tension: 0.1, pointRadius: 0 }}
                ]
            }},
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                scales: {{
                    y: {{
                        min: 0,
                        max: 100,
                        title: {{
                            display: true,
                            text: 'Percentage (%)'
                        }}
                    }},
                    x: {{
                        title: {{
                            display: true,
                            text: 'Cycle / Position (bp)'
                        }}
                    }}
                }}
            }}
        }});

        // 3. Per-Read Quality Distribution (Pie/Doughnut)
        const qBins = qcData.quality_bins;
        const qLabels = Object.keys(qBins);
        const qValues = Object.values(qBins);
        
        new Chart(document.getElementById('qualDistributionChart'), {{
            type: 'doughnut',
            data: {{
                labels: qLabels,
                datasets: [{{
                    data: qValues,
                    backgroundColor: [
                        '#ef4444', // Q < 10
                        '#f59e0b', // Q 10-19
                        '#3b82f6', // Q 20-29
                        '#10b981'  // Q >= 30
                    ],
                    borderWidth: 0
                }}]
            }},
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                plugins: {{
                    legend: {{
                        position: 'right',
                        labels: {{
                            boxWidth: 15,
                            padding: 15
                        }}
                    }}
                }},
                cutout: '65%'
            }}
        }});

        // 4. Length Distribution Chart
        const lengthLabels = Object.keys(qcData.length_distribution);
        const lengthCounts = Object.values(qcData.length_distribution);
        
        new Chart(document.getElementById('lengthChart'), {{
            type: 'bar',
            data: {{
                labels: lengthLabels,
                datasets: [{{
                    label: 'Count',
                    data: lengthCounts,
                    backgroundColor: 'rgba(6, 182, 212, 0.4)',
                    borderColor: '#06b6d4',
                    borderWidth: 1
                }}]
            }},
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                scales: {{
                    y: {{
                        beginAtZero: true,
                        title: {{
                            display: true,
                            text: 'Read Count'
                        }}
                    }},
                    x: {{
                        title: {{
                            display: true,
                            text: 'Sequence Length (bp)'
                        }}
                    }}
                }},
                plugins: {{
                    legend: {{ display: false }}
                }}
            }}
        }});
    </script>
</body>
</html>
"""

    with open(output_path, 'w') as f:
        f.write(html_content)
    print(f"Gorgeous report successfully generated at: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="FASTQ Quality Control Analyzer")
    parser.add_argument("fastq", nargs="?", help="Path to input FASTQ file (supports .fastq, .fq, or .gz versions)")
    parser.add_argument("-o", "--output", default="qc_report.html", help="Path to write the HTML report")
    parser.add_argument("--generate-sample", action="store_true", help="Generate a dummy FASTQ file for testing and exit")
    parser.add_argument("--sample-reads", type=int, default=1000, help="Number of reads to generate for sample file")
    
    args = parser.parse_args()
    
    if args.generate_sample:
        filename = args.fastq if args.fastq else "sample.fastq"
        generate_sample_fastq(filename, args.sample_reads)
        return
        
    # If no file is specified, generate a dummy sample file automatically
    fastq_file = args.fastq
    if not fastq_file:
        fastq_file = "sample.fastq"
        if not os.path.exists(fastq_file):
            print(f"No input file specified. Automatically generating sample file: {fastq_file}")
            generate_sample_fastq(fastq_file, args.sample_reads)
            
    if not os.path.exists(fastq_file):
        print(f"Error: File '{fastq_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
        
    stats = run_qc(fastq_file)
    generate_html_report(stats, args.output)

if __name__ == "__main__":
    main()
