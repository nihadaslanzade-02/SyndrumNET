"""
Automated downloaders for all required data sources.

Fetches molecular interaction networks, disease/drug expression data,
and disease susceptibility genes from public databases.
"""

import gzip
import logging
import shutil
import time
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import urljoin

import requests
from tqdm import tqdm

logger = logging.getLogger(__name__)


class DataDownloader:
    """
    Centralized data downloader for SyndrumNET pipeline.
    
    Downloads and caches data from multiple sources:
    - Molecular interactions: HuRI, CORUM, PhosphoSitePlus, etc.
    - Disease expression: CREEDS
    - Drug expression: LINCS L1000
    - Disease genes: OMIM, ClinVar, DisGeNET, etc.
    
    Parameters
    ----------
    data_dir : Path
        Root directory for downloaded data.
    timeout : int
        Request timeout in seconds.
    retry_attempts : int
        Number of retry attempts for failed downloads.
        
    Examples
    --------
    >>> downloader = DataDownloader(Path("data/raw"))
    >>> downloader.download_all()
    >>> files = downloader.get_file_paths()
    """
    
    # Data source URLs (as of January 2025)
    URLS = {
        # Molecular interactions
        "huri": "http://www.interactome-atlas.org/data/HuRI.tsv",
        "corum": "https://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip",
        "phosphositeplus": "https://www.phosphosite.org/downloads/Kinase_Substrate_Dataset.gz",
        "kegg_rpair": "https://rest.kegg.jp/get/rpair",
        "signalink": "http://signalink.org/download/edge_data",
        "innatedb": "https://www.innatedb.com/download/interactions/innatedb_ppi.mitab.gz",
        "instruct": "http://instruct.yulab.org/download/interactome.txt",
        
        # Disease expression (CREEDS)
        "creeds_disease": "https://maayanlab.cloud/Harmonizome/api/1.0/dataset/CREEDS+Disease+Signatures/download",
        
        # Drug expression (LINCS L1000)
        "lincs_l1000": "https://maayanlab.cloud/L1000CDS2/download/data/l1000_cpd_signatures.txt.gz",
        "lincs_metadata": "https://maayanlab.cloud/L1000CDS2/download/data/l1000_metadata.txt",
        
        # Disease genes
        "omim": "https://data.omim.org/downloads/[API_KEY]/genemap2.txt",  # Requires registration
        "clinvar": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
        "disgenet": "https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz",
        
        # ID mapping resources
        "hgnc": "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt",
        "uniprot_idmapping": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz",
    }
    
    def __init__(
        self,
        data_dir: Path,
        timeout: int = 300,
        retry_attempts: int = 3,
    ) -> None:
        """Initialize data downloader."""
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout
        self.retry_attempts = retry_attempts
        
        # Track versions
        self.version_file = self.data_dir / "VERSIONS.txt"
    
    def download_file(
        self,
        url: str,
        output_path: Path,
        description: Optional[str] = None,
    ) -> bool:
        """
        Download a single file with retry logic and progress bar.
        
        Parameters
        ----------
        url : str
            Source URL.
        output_path : Path
            Local destination path.
        description : str, optional
            Description for progress bar.
            
        Returns
        -------
        bool
            True if successful, False otherwise.
        """
        if output_path.exists():
            logger.info(f"File already exists: {output_path.name}")
            return True
        
        for attempt in range(self.retry_attempts):
            try:
                logger.info(f"Downloading {description or url} (attempt {attempt + 1})")
                
                response = requests.get(url, stream=True, timeout=self.timeout)
                response.raise_for_status()
                
                # Get file size if available
                total_size = int(response.headers.get('content-length', 0))
                
                # Download with progress bar
                output_path.parent.mkdir(parents=True, exist_ok=True)
                
                with open(output_path, 'wb') as f:
                    with tqdm(
                        total=total_size,
                        unit='B',
                        unit_scale=True,
                        desc=output_path.name,
                    ) as pbar:
                        for chunk in response.iter_content(chunk_size=8192):
                            f.write(chunk)
                            pbar.update(len(chunk))
                
                logger.info(f"Successfully downloaded: {output_path.name}")
                return True
                
            except Exception as e:
                logger.warning(f"Download attempt {attempt + 1} failed: {e}")
                if attempt < self.retry_attempts - 1:
                    time.sleep(2 ** attempt)  # Exponential backoff
                else:
                    logger.error(f"Failed to download {url} after {self.retry_attempts} attempts")
                    return False
        
        return False
    
    def download_huri(self) -> Path:
        """Download HuRI human interactome."""
        output = self.data_dir / "huri.tsv"
        self.download_file(self.URLS["huri"], output, "HuRI")
        return output
    
    def download_corum(self) -> Path:
        """Download CORUM protein complexes."""
        output = self.data_dir / "corum.zip"
        self.download_file(self.URLS["corum"], output, "CORUM")
        
        # Extract if downloaded
        if output.exists():
            import zipfile
            with zipfile.ZipFile(output, 'r') as zf:
                zf.extractall(self.data_dir)
            logger.info("Extracted CORUM data")
        
        return self.data_dir / "allComplexes.txt"
    
    def download_phosphositeplus(self) -> Path:
        """
        Download PhosphoSitePlus kinase-substrate data.
        
        Note: Requires registration at phosphosite.org
        Downloads may require manual authentication.
        """
        output = self.data_dir / "phosphositeplus.txt.gz"
        
        logger.warning(
            "PhosphoSitePlus requires registration. "
            "If download fails, manually download from phosphosite.org "
            "and place in data/raw/"
        )
        
        self.download_file(self.URLS["phosphositeplus"], output, "PhosphoSitePlus")
        
        # Decompress
        if output.exists():
            decompressed = output.with_suffix('')
            with gzip.open(output, 'rb') as f_in:
                with open(decompressed, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            return decompressed
        
        return output
    
    def download_creeds(self) -> Path:
        """Download CREEDS disease signatures."""
        output = self.data_dir / "creeds_disease.txt"
        self.download_file(self.URLS["creeds_disease"], output, "CREEDS")
        return output
    
    def download_lincs(self) -> Dict[str, Path]:
        """Download LINCS L1000 drug signatures and metadata."""
        sig_file = self.data_dir / "lincs_signatures.txt.gz"
        meta_file = self.data_dir / "lincs_metadata.txt"
        
        self.download_file(self.URLS["lincs_l1000"], sig_file, "LINCS L1000 signatures")
        self.download_file(self.URLS["lincs_metadata"], meta_file, "LINCS metadata")
        
        # Decompress signatures
        if sig_file.exists():
            decompressed = sig_file.with_suffix('')
            with gzip.open(sig_file, 'rb') as f_in:
                with open(decompressed, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            return {
                "signatures": decompressed,
                "metadata": meta_file,
            }
        
        return {"signatures": sig_file, "metadata": meta_file}
    
    def download_disease_genes(self) -> Dict[str, Path]:
        """Download disease-gene associations from multiple sources."""
        files = {}
        
        # ClinVar
        clinvar = self.data_dir / "clinvar.txt.gz"
        if self.download_file(self.URLS["clinvar"], clinvar, "ClinVar"):
            files["clinvar"] = clinvar
        
        # DisGeNET
        disgenet = self.data_dir / "disgenet.tsv.gz"
        if self.download_file(self.URLS["disgenet"], disgenet, "DisGeNET"):
            files["disgenet"] = disgenet
        
        # OMIM requires API key
        logger.warning(
            "OMIM download requires API key from omim.org. "
            "Set API key in config or download manually."
        )
        
        return files
    
    def download_id_mapping(self) -> Dict[str, Path]:
        """Download gene/protein ID mapping resources."""
        files = {}
        
        # HGNC
        hgnc = self.data_dir / "hgnc.txt"
        if self.download_file(self.URLS["hgnc"], hgnc, "HGNC"):
            files["hgnc"] = hgnc
        
        # UniProt ID mapping
        uniprot = self.data_dir / "uniprot_idmapping.dat.gz"
        if self.download_file(self.URLS["uniprot_idmapping"], uniprot, "UniProt"):
            files["uniprot"] = uniprot
        
        return files
    
    def download_all(self) -> Dict[str, Path]:
        """
        Download all required data sources.
        
        Returns
        -------
        dict
            Mapping of source names to file paths.
        """
        logger.info("Starting download of all data sources")
        
        files = {}
        
        # Molecular interactions
        logger.info("Downloading molecular interaction data...")
        files["huri"] = self.download_huri()
        files["corum"] = self.download_corum()
        files["phosphositeplus"] = self.download_phosphositeplus()
        
        # Disease and drug expression
        logger.info("Downloading expression data...")
        files["creeds"] = self.download_creeds()
        files.update(self.download_lincs())
        
        # Disease genes
        logger.info("Downloading disease-gene associations...")
        files.update(self.download_disease_genes())
        
        # ID mapping
        logger.info("Downloading ID mapping resources...")
        files.update(self.download_id_mapping())
        
        # Save version info
        self._save_versions(files)
        
        logger.info(f"Downloaded {len(files)} data files")
        return files
    
    def _save_versions(self, files: Dict[str, Path]) -> None:
        """Save version/date information for downloaded files."""
        from datetime import datetime
        
        with open(self.version_file, 'w') as f:
            f.write(f"# Data versions - Downloaded on {datetime.now()}\n")
            f.write(f"# SyndrumNET Reproduction Pipeline\n\n")
            
            for source, path in files.items():
                if path.exists():
                    size = path.stat().st_size / (1024 * 1024)  # MB
                    f.write(f"{source}: {path.name} ({size:.2f} MB)\n")
        
        logger.info(f"Saved version info to {self.version_file}")
    
    def get_file_paths(self) -> Dict[str, Path]:
        """Get paths to all expected data files."""
        return {
            "huri": self.data_dir / "huri.tsv",
            "corum": self.data_dir / "allComplexes.txt",
            "creeds": self.data_dir / "creeds_disease.txt",
            "lincs_sigs": self.data_dir / "lincs_signatures.txt",
            "lincs_meta": self.data_dir / "lincs_metadata.txt",
            "hgnc": self.data_dir / "hgnc.txt",
        }