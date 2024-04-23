use polars::prelude::*;
use clap::{Command, Arg};
use std::path::PathBuf;
use std::io::prelude::*;
use std::io::{BufWriter};
use std::fs::File;


fn cli() -> Command<'static> {
    Command::new("downsample")
             .args(&[
                 Arg::new("input")
                      .long("input")
                      .short('i')
                      .takes_value(true)
                      .required(true)
                      .value_name("FILE"),
                 Arg::new("output")
                      .long("output")
                      .short('o')
                      .takes_value(true)
                      .required(true)
                      .value_name("FILE"),
                 Arg::new("fraction")
                      .long("frac")
                      .short('f')
                      .takes_value(true)
                      .required(true)
                      .value_name("F32"),
                 Arg::new("seed")
                      .long("seed")
                      .short('s')
                      .takes_value(true)
                      .value_name("U64")
             ])
}

fn main() -> Result<()> {
     let matches = cli().get_matches();

     let input = matches.value_of("input")
                         .unwrap();
     let output = matches.value_of("output")
                         .unwrap();
     let seed = matches.value_of("seed")
                         .map(|v| 
                              v.parse::<u64>()
                              .expect("seed was not a valid u64")
                         );
     let frac = matches.value_of("fraction")
                         .map(|v| 
                              v.parse::<f64>()
                              .expect("frac was not an f64"))
                         .unwrap();

     let mut bedpe = CsvReader::from_path(PathBuf::from(input.to_string()))?
                              .infer_schema(Some(5))
                              .with_delimiter(b'\t')
                              .has_header(false)
                              .with_encoding(CsvEncoding::LossyUtf8)
                              .low_memory(true)
                              .finish()
                              .expect("failed to read csv");

     let keep_chroms = Series::new("keep_chroms", &[
          "chr1",
          "chr2",
          "chr3",
          "chr4",
          "chr5",
          "chr6",
          "chr7",
          "chrX",
          "chr8",
          "chr9",
          "chr11",
          "chr10",
          "chr12",
          "chr13",
          "chr14",
          "chr15",
          "chr16",
          "chr17",
          "chr18",
          "chr20",
          "chr19",
          "chrY",
          "chr22",
          "chr21"
      ]);

      bedpe.set_column_names(&["chrom1", "start1", "end1", 
                               "chrom2", "start2", "end2",
                               "name", "score", "strand1", "strand2"])?;

     let bedpe = bedpe.filter(
                         &bedpe["chrom1"].is_in(&keep_chroms).expect("Error masking good chromosomes")
               ).expect("Error filtering good chromosomes");

     let mut bedpe = bedpe.sample_frac(frac,
                                   false,
                                   seed)?;

     
     bedpe.sort_in_place(&["chrom1", "start1"],
                         vec![false, false])?;

     let buf_out = BufWriter::with_capacity(100 * 1024_usize.pow(2),
                                             File::create(PathBuf::from(output.to_string())).unwrap());

     
     let _result = CsvWriter::new(buf_out)
                              .with_delimiter(b'\t')
                              .has_header(false)
                              .finish(&mut bedpe)
                              .expect("failed to write csv");
     Ok(())
}
