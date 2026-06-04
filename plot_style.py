"""SciencePlots nature style for all project figures."""
import matplotlib.pyplot as plt


def apply():
    import scienceplots  # noqa: F401 — registers styles with matplotlib

    plt.style.use(["science", "nature", "no-latex"])
