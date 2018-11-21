CREATE TABLE gene
  (
    entrez_id integer,
    symbol varchar,
    ensembl_gene_id varchar
  );

CREATE TABLE variant
  (
    name varchar,
    biomarker_type varchar,
    so_id varchar,
    so_hierarchy varchar[]
  );
