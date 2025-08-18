from pathlib import Path
import time
import subprocess
from tools.functions import setup_logger
logger = setup_logger()


class PerturbationManager:
    """
    Thread-safe perturbation manager using file-based locking.
    Ensures only one process perturbs a dataset at a time.
    """
    
    def __init__(self, base_dir: str = "output"):
        self.base_dir = Path(base_dir)
        self.lock_dir = self.base_dir / "locks"
        self.lock_dir.mkdir(parents=True, exist_ok=True)

    def get_perturbation_key(self, dataset: str, change_percent: float, iteration: int) -> str:
        """Create a unique key for each perturbation configuration."""
        return f"{dataset}_{change_percent:.4f}_iter{iteration}"
    
    def get_lock_file(self, key: str) -> Path:
        """Get the lock file path for a perturbation key."""
        return self.lock_dir / f"{key}.lock"
    
    def get_done_file(self, key: str) -> Path:
        """Get the completion marker file path."""
        return self.lock_dir / f"{key}.done"
    
    def is_perturbation_done(self, dataset: str, change_percent: float, iteration: int) -> bool:
        """Check if perturbation is already completed."""
        key = self.get_perturbation_key(dataset, change_percent, iteration)
        done_file = self.get_done_file(key)
        return done_file.exists()
    
    def run_perturbation_safe(self, dataset: str, change_percent: float, iteration: int,  
                             script_path: str, timeout: int = 300) -> bool:
        """
        Run perturbation in a thread-safe manner using file locking.
        """
        key = self.get_perturbation_key(dataset, change_percent, iteration)
        lock_file = self.get_lock_file(key)
        done_file = self.get_done_file(key)
        
        # Check if already done
        if done_file.exists():
            logger.info(f"Perturbation already completed: {key}")
            return True
        
        # Try to acquire lock
        max_wait_time = 600  # 10 minutes max wait
        wait_time = 0
        check_interval = 5  # Check every 5 seconds
        
        while wait_time < max_wait_time:
            try:
                # Try to create lock file atomically
                lock_file.touch(exist_ok=False)
                logger.info(f"Acquired lock for perturbation: {key}")
                break
            except FileExistsError:
                # Another process has the lock, wait
                logger.info(f"Waiting for perturbation lock: {key} (waited {wait_time}s)")
                time.sleep(check_interval)
                wait_time += check_interval
                
                # Check if perturbation completed while waiting
                if done_file.exists():
                    logger.info(f"Perturbation completed by another process: {key}")
                    return True
        else:
            # Timeout waiting for lock
            logger.error(f"Timeout waiting for perturbation lock: {key}")
            return False
        
        try:
            # Double-check if done (race condition protection)
            if done_file.exists():
                logger.info(f"Perturbation already completed by another process: {key}")
                return True
            
            # Run the actual perturbation
            logger.info(f"Running perturbation: {dataset} ({change_percent:.1%}), iteration={iteration}")         

            # Calculate seed based on iteration
            seed = 42 * iteration
            
            cmd = [
                "Rscript", script_path,
                "-d", dataset,
                "-p", str(change_percent),
                "-s", str(seed)
            ]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout
            )
            
            if result.returncode != 0:
                logger.error(f"Perturbation failed for {key}: {result.stderr}")
                return False
            
            # Mark as completed
            done_file.touch()
            logger.info(f"Perturbation completed successfully: {key}")
            return True
            
        except subprocess.TimeoutExpired:
            logger.error(f"Perturbation timeout for {key}")
            return False
        except Exception as e:
            logger.error(f"Perturbation error for {key}: {e}")
            return False
        finally:
            # Always release lock
            try:
                if lock_file.exists():
                    lock_file.unlink()
                    logger.debug(f"Released lock for: {key}")
            except Exception as e:
                logger.warning(f"Failed to release lock for {key}: {e}")
    
    def cleanup_locks(self, dataset: str = None):
        """Clean up lock files, optionally for a specific dataset."""
        try:
            if dataset:
                # Clean locks for specific dataset
                pattern = f"{dataset}_*.lock"
                lock_files = list(self.lock_dir.glob(pattern))
                done_pattern = f"{dataset}_*.done"
                done_files = list(self.lock_dir.glob(done_pattern))
                files_to_remove = lock_files + done_files
            else:
                # Clean all locks
                files_to_remove = list(self.lock_dir.glob("*.lock")) + list(self.lock_dir.glob("*.done"))
            
            for file_path in files_to_remove:
                try:
                    file_path.unlink()
                    logger.debug(f"Removed lock file: {file_path}")
                except Exception as e:
                    logger.warning(f"Failed to remove {file_path}: {e}")
                    
        except Exception as e:
            logger.error(f"Error during lock cleanup: {e}")