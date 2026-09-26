__version__ = "0.4.0"

from .engine import (
    TideConfig, TideFit, TideResult, fit_historical, test_treatment, get_envelope,
    mann_kendall_test, sens_slope, to_day_of_year, assign_bins, n_bins_for,
    holm_bonferroni, bin_significant, critical_value, significance_band, SustainedResult,
    DepartureRecovery, estimate_departure_recovery,
)
from .controllers import (
    SequentialEvent, SequentialResultRow, sequential_test,
    CumulativeWindow, CumulativeResultRow, cumulative_test,
    RecoveryResult, first_recovery,
)
from .power import estimate_power, estimate_power_sequential, sensitivity_grid, simulate_synthetic_curve
from .monitoring import MonitoringSeries, MonitoringCheck
from .plotting import plot_envelope, plot_sequential, plot_cumulative, plot_monitoring
from .report import summarize

__all__ = [
    "__version__",
    "TideConfig", "TideFit", "TideResult", "fit_historical", "test_treatment", "get_envelope",
    "mann_kendall_test", "sens_slope", "to_day_of_year", "assign_bins", "n_bins_for",
    "holm_bonferroni", "bin_significant", "critical_value", "significance_band", "SustainedResult",
    "DepartureRecovery", "estimate_departure_recovery",
    "SequentialEvent", "SequentialResultRow", "sequential_test",
    "CumulativeWindow", "CumulativeResultRow", "cumulative_test",
    "RecoveryResult", "first_recovery",
    "estimate_power", "estimate_power_sequential", "sensitivity_grid", "simulate_synthetic_curve",
    "MonitoringSeries", "MonitoringCheck",
    "plot_envelope", "plot_sequential", "plot_cumulative", "plot_monitoring",
    "summarize",
]
