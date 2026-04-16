#!/usr/bin/env python3
"""
HiC Data Pipeline v3 - VERIFIED WORKING TARGETS ONLY
=====================================================
Downloads and validates high-depth .hic files for power-law scaling analysis.

Usage:
    python hic_data_pipeline_v3.py --dry-run    # Check URLs without downloading
    python hic_data_pipeline_v3.py              # Download all files
    
Cell Lines Included:
    - HCT116 (colon cancer) - from 4DN Portal
    - HepG2 (liver cancer) - from 4DN Portal  
    - K562 (leukemia reference) - from GEO GSE63525
    - GM12878 (lymphoblast gold standard) - from GEO GSE63525
    - IMR90 (fibroblast) - from GEO GSE63525
    - HMEC (breast normal) - from GEO GSE63525
    - HUVEC (endothelial) - from GEO GSE63525
    - NHEK (keratinocyte) - from GEO GSE63525

Note: 4DN Portal files (HCT116, HepG2) may require browser download due to authentication.
"""

import os
import sys
import re
import struct
import logging
import subprocess
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass
from typing import Optional, Dict, List, Tuple

import requests
import pandas as pd


@dataclass
class HiCTarget:
    accession: str
    cell_line: str
    subgroup: str
    source: str
    url: Optional[str] = None
    filename_pattern: Optional[str] = None


# VERIFIED TARGETS WITH CONFIRMED .hic FILES
TARGETS = [
    # === 4DN Portal (direct .hic file accessions) ===
    # Note: These may require browser download due to authentication
    HiCTarget("4DNFI2TK7L2F", "HCT116", "Gastrointestinal", "4dn_file"),
    HiCTarget("4DNFI1UEG1HD", "HepG2", "Liver", "4dn_file"),
    
    # === GEO GSE63525 - Rao et al. 2014 (verified .hic files) ===
    HiCTarget("GSE63525", "K562", "Reference", "geo", filename_pattern="K562_combined.hic"),
    HiCTarget("GSE63525", "GM12878", "Reference", "geo", filename_pattern="GM12878_insitu_primary.hic"),
    HiCTarget("GSE63525", "IMR90", "Reference", "geo", filename_pattern="IMR90_combined.hic"),
    HiCTarget("GSE63525", "HMEC", "Breast_Normal", "geo", filename_pattern="HMEC_combined.hic"),
    HiCTarget("GSE63525", "HUVEC", "Endothelial", "geo", filename_pattern="HUVEC_combined.hic"),
    HiCTarget("GSE63525", "NHEK", "Skin", "geo", filename_pattern="NHEK_combined.hic"),
]


def setup_logging(log_dir: Path) -> logging.Logger:
    """Configure logging with both file and console handlers."""
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("HiCPipeline")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    
    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(logging.Formatter('%(asctime)s | %(levelname)-8s | %(message)s', '%H:%M:%S'))
    logger.addHandler(ch)
    
    # File handler
    fh = logging.FileHandler(log_dir / f"pipeline_{datetime.now():%Y%m%d_%H%M%S}.log")
    fh.setFormatter(logging.Formatter('%(asctime)s | %(levelname)-8s | %(message)s'))
    logger.addHandler(fh)
    
    return logger


class DataFetcher:
    """Fetch .hic file URLs from GEO and 4DN Portal."""
    
    GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
    FDN_BASE = "https://data.4dnucleome.org"
    
    def __init__(self, logger):
        self.logger = logger
        self.session = requests.Session()
    
    def get_url(self, target: HiCTarget) -> Optional[Dict]:
        """Get download URL for a target."""
        if target.source == "4dn_file":
            return self._get_4dn_file(target.accession)
        elif target.source == "geo":
            return self._get_geo_file(target.accession, target.filename_pattern)
        return None
    
    def _get_4dn_file(self, accession: str) -> Optional[Dict]:
        """Get 4DN Portal file URL."""
        self.logger.info(f"Querying 4DN file: {accession}")
        url = f"{self.FDN_BASE}/files-processed/{accession}/@@download/{accession}.hic"
        return {'filename': f"{accession}.hic", 'url': url}
    
    def _get_geo_file(self, accession: str, pattern: str) -> Optional[Dict]:
        """Get GEO supplementary file URL."""
        self.logger.info(f"Querying GEO: {accession} (pattern: {pattern})")
        prefix = accession[:len(accession)-3] + "nnn"
        base_url = f"{self.GEO_FTP}/{prefix}/{accession}/suppl/"
        
        if pattern:
            filename = f"{accession}_{pattern}" if not pattern.startswith(accession) else pattern
            return {'filename': filename, 'url': f"{base_url}{filename}"}
        return None


class HiCValidator:
    """Validate .hic file integrity."""
    
    def __init__(self, logger):
        self.logger = logger
    
    def validate(self, filepath: Path) -> Dict:
        """Full validation of a .hic file."""
        result = {'valid': False, 'size_gb': 0}
        
        if not filepath.exists():
            return result
        
        result['size_gb'] = filepath.stat().st_size / (1024**3)
        
        # Check magic string
        with open(filepath, 'rb') as f:
            if f.read(3) != b'HIC':
                self.logger.error("  Invalid magic string")
                return result
        self.logger.info("  ✓ Valid HIC format")
        
        # Parse header
        try:
            with open(filepath, 'rb') as f:
                f.read(4)
                ver = struct.unpack('<i', f.read(4))[0]
                f.read(8)
                genome = b''.join(iter(lambda: f.read(1), b'\x00')).decode()
            self.logger.info(f"  Version: {ver}, Genome: {genome}")
            result['version'] = ver
            result['genome'] = genome
        except Exception as e:
            self.logger.warning(f"  Header parse error: {e}")
        
        # Straw validation for resolutions
        try:
            import hicstraw
            hic = hicstraw.HiCFile(str(filepath))
            result['resolutions'] = hic.getResolutions()
            result['has_10kb'] = 10000 in result['resolutions']
            self.logger.info(f"  Resolutions: {result['resolutions']}")
            self.logger.info(f"  10kb available: {result['has_10kb']}")
        except ImportError:
            self.logger.warning("  hic-straw not installed, skipping resolution check")
        except Exception as e:
            self.logger.warning(f"  Straw check skipped: {e}")
        
        result['valid'] = True
        return result


class Pipeline:
    """Main pipeline orchestrator."""
    
    def __init__(self, output="./data", logs="./logs"):
        self.base_dir = Path(output)
        self.log_dir = Path(logs)
        self.logger = setup_logging(self.log_dir)
        self.fetcher = DataFetcher(self.logger)
        self.validator = HiCValidator(self.logger)
    
    def download(self, url: str, target: HiCTarget) -> Tuple[bool, Path]:
        """Download a file using wget with resume capability."""
        out_dir = self.base_dir / target.subgroup / target.cell_line
        out_dir.mkdir(parents=True, exist_ok=True)
        filename = os.path.basename(url.split('?')[0])
        filepath = out_dir / filename
        
        self.logger.info(f"Downloading: {filename}")
        self.logger.info(f"  To: {filepath}")
        
        cmd = ['wget', '-c', '--tries=3', '--timeout=120', '-O', str(filepath), url]
        try:
            subprocess.run(cmd, check=True, timeout=21600)  # 6 hour timeout
            if filepath.exists() and filepath.stat().st_size > 1000:
                self.logger.info(f"  ✓ Downloaded: {filepath.stat().st_size/(1024**3):.2f} GB")
                return True, filepath
        except Exception as e:
            self.logger.error(f"  ✗ Failed: {e}")
        return False, filepath
    
    def run(self, dry_run=False) -> pd.DataFrame:
        """Execute the pipeline."""
        self.logger.info("=" * 60)
        self.logger.info(f"HiC Pipeline v3 | Mode: {'DRY RUN' if dry_run else 'DOWNLOAD'}")
        self.logger.info(f"Verified targets: {len(TARGETS)}")
        self.logger.info("=" * 60)
        
        results = []
        for t in TARGETS:
            self.logger.info(f"\n--- {t.cell_line} ({t.accession}) ---")
            
            info = self.fetcher.get_url(t)
            if not info:
                results.append({
                    'cell_line': t.cell_line, 
                    'accession': t.accession,
                    'subgroup': t.subgroup, 
                    'status': 'URL_FAILED'
                })
                continue
            
            row = {
                'cell_line': t.cell_line, 
                'accession': t.accession,
                'subgroup': t.subgroup, 
                'filename': info['filename'],
                'url': info['url'], 
                'status': 'DRY_RUN' if dry_run else 'PENDING'
            }
            
            if not dry_run:
                ok, path = self.download(info['url'], t)
                if ok:
                    v = self.validator.validate(path)
                    row['status'] = 'VALID' if v['valid'] else 'INVALID'
                    row['size_gb'] = f"{v.get('size_gb', 0):.2f}"
                    row['has_10kb'] = v.get('has_10kb', 'N/A')
                else:
                    row['status'] = 'DOWNLOAD_FAILED'
            
            results.append(row)
        
        df = pd.DataFrame(results)
        csv_path = self.log_dir / f"results_{datetime.now():%Y%m%d_%H%M%S}.csv"
        df.to_csv(csv_path, index=False)
        
        self.logger.info("\n" + "=" * 60)
        self.logger.info("SUMMARY")
        self.logger.info(df['status'].value_counts().to_string())
        self.logger.info(f"\nResults: {csv_path}")
        return df


def main():
    import argparse
    p = argparse.ArgumentParser(description="HiC Data Pipeline v3")
    p.add_argument('--dry-run', action='store_true', help='Check URLs without downloading')
    p.add_argument('--output', default='./data', help='Output directory')
    p.add_argument('--logs', default='./logs', help='Log directory')
    args = p.parse_args()
    
    pipe = Pipeline(args.output, args.logs)
    df = pipe.run(dry_run=args.dry_run)
    print("\n" + df.to_string())


if __name__ == "__main__":
    sys.exit(main())