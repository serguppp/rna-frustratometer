import csv
import math
from itertools import product
from statistics import fmean, pstdev

import pytest

from rna_frustratometer.frustration import (
    COMBINATIONS,
    BaseFrequencies,
    energy_statistics,
    save_frustration_scores,
)
from rna_frustratometer.pipeline import ContactResult


def test_2koc_exact_matrix_and_uniform_baseline():
    background = BaseFrequencies.from_sequence("GGCACUUCGGUGCC")
    assert background.counts == (1, 5, 5, 3)
    assert background.weights() == pytest.approx(
        [n / 196 for n in (1, 5, 5, 3, 5, 25, 25, 15, 5, 25, 25, 15, 3, 15, 15, 9)]
    )
    assert sum(background.weights()) == pytest.approx(1)
    assert background.weights(mode="uniform") == (1 / 16,) * 16


def test_weighted_moments_match_independent_explicit_sum():
    weights = BaseFrequencies.from_sequence("ACCCCCGGUUU").weights()
    energies = [float(k * k - 3 * k) for k in range(16)]
    mean = sum(w * e for w, e in zip(weights, energies, strict=True))
    sd = math.sqrt(sum(w * (e - mean) ** 2 for w, e in zip(weights, energies, strict=True)))
    result = energy_statistics(-5, energies, weights)
    assert result.mean == pytest.approx(mean)
    assert result.sd == pytest.approx(sd)
    assert result.frustration == pytest.approx((mean + 5) / sd)


def test_uniform_matches_previous_population_statistics_and_shift_invariance():
    energies = [float(k) for k in range(16)]
    result = energy_statistics(2, energies)
    assert result.mean == fmean(energies)
    assert result.sd == pstdev(energies)
    assert result.frustration == (fmean(energies) - 2) / pstdev(energies)
    shifted = energy_statistics(10002, [10000 + e for e in energies])
    assert shifted.frustration == result.frustration


def test_zero_probability_states_and_optional_pseudocount():
    background = BaseFrequencies.from_sequence("AAAA")
    assert background.weights() == (1.0,) + (0.0,) * 15
    assert energy_statistics(3, list(range(16)), background.weights()).frustration is None
    regularized = BaseFrequencies.from_sequence("AAAA", pseudocount=1)
    assert regularized.probabilities == (5 / 8, 1 / 8, 1 / 8, 1 / 8)
    assert min(regularized.weights()) > 0


@pytest.mark.parametrize("bad", ["", "ACNT", "ACT", "AX"])
def test_ambiguous_or_empty_sequence_is_rejected(bad):
    with pytest.raises(ValueError):
        BaseFrequencies.from_sequence(bad)


def test_invalid_weights_and_numerical_resolution_are_rejected_or_undefined():
    with pytest.raises(ValueError):
        energy_statistics(0, list(range(16)), [1] * 16)
    with pytest.raises(ValueError):
        energy_statistics(0, [float("nan")] * 16)
    assert energy_statistics(4, [4] * 16).frustration is None
    assert (
        energy_statistics(4, [4 + k * 1e-8 for k in range(16)], minimum_sd=2e-6).frustration is None
    )


def test_variant_order_does_not_change_weights_or_result():
    background = BaseFrequencies.from_sequence("ACCCGGUUUU")
    energies = tuple(float(i) for i in range(16))
    first = energy_statistics(0, energies, background.weights())
    second = energy_statistics(0, energies[::-1], background.weights(COMBINATIONS[::-1]))
    assert first.frustration == pytest.approx(second.frustration)


def test_contact_result_exposes_configured_and_uniform_scores():
    energies = tuple(float(i) for i in range(16))
    weights = BaseFrequencies.from_sequence("ACCCGGUUUU").weights()
    contact = ContactResult("A", "1", "A", "2", 0, energies, weights)
    assert contact.frustration == energy_statistics(0, energies, weights).frustration
    assert contact.baseline_frustration == energy_statistics(0, energies).frustration


def test_saved_scores_include_both_strategies_and_alias_defaults_to_composition(tmp_path):
    source = tmp_path / "energies.csv"
    with source.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "chain_i",
                "residue_i",
                "chain_j",
                "residue_j",
                "base_i",
                "base_j",
                "native_E_pair",
                "E_pair",
                "native_E_total",
                "E_total",
            ]
        )
        for index, (a, b) in enumerate(product("ACGU", repeat=2)):
            writer.writerow(["A", 1, "A", 2, a, b, 0, index, 20, 20 + 2 * index])
    output = save_frustration_scores(
        source, tmp_path / "scores", BaseFrequencies.from_sequence("ACCCGGUUUU")
    )
    with output.open(newline="", encoding="utf-8") as stream:
        row = next(csv.DictReader(stream))
    assert row["F_pair"] == row["F_pair_composition"]
    assert row["F_total"] == row["F_total_composition"]
    assert float(row["F_pair_uniform"]) == pytest.approx(fmean(range(16)) / pstdev(range(16)))
    assert row["F_pair_composition"] != row["F_pair_uniform"]
    assert (tmp_path / "scores/decoy_weights.json").is_file()


def test_pdb_frequencies_count_residues_once_not_atoms_or_nmr_models(tmp_path):
    import gemmi

    structure = gemmi.Structure()
    for model_number in (1, 2):
        model = gemmi.Model(str(model_number))
        chain = gemmi.Chain("A")
        for number, base in enumerate("ACGU", 1):
            residue = gemmi.Residue()
            residue.name = base
            residue.seqid = gemmi.SeqId(number, " ")
            for name in ("P", "C4'"):
                atom = gemmi.Atom()
                atom.name = name
                residue.add_atom(atom)
            chain.add_residue(residue)
        model.add_chain(chain)
        structure.add_model(model)
    path = tmp_path / "two_models.pdb"
    structure.write_pdb(str(path))
    background = BaseFrequencies.from_pdb(path)
    assert background.counts == (1, 1, 1, 1)
    assert background.weights() == (1 / 16,) * 16
