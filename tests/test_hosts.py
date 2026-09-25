"""Host names from PLSDB, PlasmidScope and GenBank/RefSeq, reduced to species and genus."""
from plasmidann import hosts


def test_a_binomial_gives_species_and_genus():
    assert hosts.normalise("Escherichia coli") == ("Escherichia coli", "Escherichia")
    # PLSDB writes underscores; a subspecies or strain is the same species.
    assert hosts.normalise("Klebsiella_pneumoniae") == ("Klebsiella pneumoniae", "Klebsiella")
    assert hosts.normalise("Salmonella enterica subsp. enterica serovar Typhimurium") == \
        ("Salmonella enterica", "Salmonella")


def test_a_genus_level_name_counts_only_toward_the_genus():
    assert hosts.normalise("Acidovorax sp. DW039") == ("", "Acidovorax")
    assert hosts.normalise("Candidatus Azobacteroides pseudotrichonymphae") == \
        ("Candidatus Azobacteroides pseudotrichonymphae", "Candidatus Azobacteroides")


def test_a_name_that_is_not_an_organism_is_no_host():
    for name in ("human gut metagenome", "uncultured bacterium", "Plasmid pDS56", "-", "",
                 "uncultured prokaryote", "bacterium"):
        assert hosts.normalise(name) == ("", ""), name


def test_the_first_source_that_names_an_organism_wins():
    assert hosts.resolve([("plsdb", ""), ("plasmidscope", "gut metagenome"),
                          ("organism", "Escherichia coli")]) == \
        ("Escherichia coli", "Escherichia", "organism")
    assert hosts.resolve([("plsdb", "Shigella_flexneri"),
                          ("organism", "Escherichia coli")])[2] == "plsdb"
    assert hosts.resolve([("plsdb", "")]) == ("", "", "")
