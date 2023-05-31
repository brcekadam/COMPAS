from compas_python_utils.cosmic_integration.totalMassEvolvedPerZ import (
    IMF, get_COMPAS_fraction
)
import numpy as np
import matplotlib.pyplot as plt


def test_imf(test_archive_dir):
    m = np.geomspace(5e-3, 200, 100)
    m_breaks = [0.01, 0.08, 0.5, 200]
    powers = [0.3, 1.3, 2.3]
    imf_values = IMF(m, *m_breaks, *powers)

    assert imf_values[0] == 0, "IMF(m<min-m) ==0"
    assert np.sum(imf_values < 0) == 0, "IMF must be > 0"

    plt.plot(m, imf_values)
    for mi in m_breaks:
        plt.axvline(mi, zorder=-10, color='gray', alpha=0.2)
    plt.xscale("log")
    plt.xlabel(r"Mass [M$_{\odot}]$")
    plt.ylabel("IMF")
    plt.savefig(f"{test_archive_dir}/IMF.png")


def test_compas_fraction():
    assert get_COMPAS_fraction(
        m1_low=10, m1_upp=100, m2_low=10, f_bin=0
    ) == 0
    assert get_COMPAS_fraction(
        m1_low=0.01, m1_upp=200, m2_low=0, f_bin=1,
        m1=0.01, m2=0.08, m3=0.5, m4=200
    ) == 1
    assert 0 < (get_COMPAS_fraction(
        m1_low=0.01, m1_upp=100, m2_low=0, f_bin=0.5,
        m1=0.01, m2=0.08, m3=0.5, m4=200
    )) < 1
