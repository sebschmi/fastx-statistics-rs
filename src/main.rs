use clap::Parser;
use indicatif::{ProgressBar, ProgressState, ProgressStyle};
use log::{info, LevelFilter};
use seq_io::fastx::Reader;
use seq_io::BaseRecord;
use simplelog::{ColorChoice, CombinedLogger, Config, TermLogger, TerminalMode};
use std::fmt::Write;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;
use std::time::{Duration, Instant};

#[derive(Parser)]
struct Cli {
    /// Fasta or fastq input file (automatically detected).
    #[clap(index = 1)]
    input: PathBuf,

    /// Filter fasta or fastq records with the given ids (pass multiple times for multiple ids).
    #[clap(long = "filter-id", value_name = "FILTER_ID")]
    filter_ids: Vec<String>,

    /// Additional percentiles for Nx metric to report (pass multiple times for multiple percentiles).
    /// N50 and N75 are always reported.
    #[clap(long = "additional-percentile", value_name = "ADDITIONAL_PERCENTILE")]
    additional_percentiles: Vec<u8>,
}

pub fn initialise_logging(log_level: LevelFilter) {
    CombinedLogger::init(vec![TermLogger::new(
        if cfg!(debug_assertions) {
            LevelFilter::Trace
        } else {
            log_level
        },
        Config::default(),
        TerminalMode::Mixed,
        ColorChoice::Auto,
    )])
    .unwrap();

    info!("Logging initialised successfully");
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();
    initialise_logging(LevelFilter::Info);

    if !cli.input.is_file() {
        return Err(format!("Not a file: {:?}", cli.input));
    }

    let input_file =
        File::open(&cli.input).map_err(|err| format!("Cannot open input file: {}", err))?;
    let input_len = input_file
        .metadata()
        .map_err(|err| format!("Cannot read file metadata: {}", err))?
        .len();
    basic_statistics(
        input_file,
        input_len,
        &cli.filter_ids,
        &cli.additional_percentiles,
    )
}

fn basic_statistics(
    input: impl Read,
    input_len: u64,
    filter_ids: &[String],
    additional_percentiles: &[u8],
) -> Result<(), String> {
    let mut fastx_reader = Reader::new(BufReader::new(input));

    let mut sequence_lengths = Vec::new();
    let mut sequence_hoco_lengths = Vec::new();
    let mut sequence_lengths_without_ns = Vec::new();
    let mut sequence_hoco_lengths_without_ns = Vec::new();

    info!("Reading fasta or fastq file...");
    let pb = ProgressBar::new(input_len);
    pb.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes} ({eta})",
        )
        .unwrap()
        .with_key("eta", |state: &ProgressState, w: &mut dyn Write| {
            write!(w, "{:.0}s", state.eta().as_secs_f64()).unwrap()
        })
        .progress_chars("#>-"),
    );
    let mut last_update = Instant::now();

    while let Some(record) = fastx_reader.next() {
        let record = record.map_err(|err| format!("Error parsing fastx: {}", err))?;
        if filter_ids.contains(
            &record
                .id()
                .map_err(|err| format!("Record id is not utf-8 encoded: {err}"))?
                .to_owned(),
        ) {
            continue;
        }

        let sequence_statistics = SequenceStatistics::new(record.seq());

        sequence_lengths.push(sequence_statistics.len);
        sequence_hoco_lengths.push(sequence_statistics.hoco_len);
        sequence_lengths_without_ns.push(sequence_statistics.len_without_ns);
        sequence_hoco_lengths_without_ns.push(sequence_statistics.hoco_len_without_ns);

        let now = Instant::now();
        if last_update + Duration::from_millis(200) <= now {
            pb.set_position(fastx_reader.position().byte());
            last_update = now;
        }
    }

    pb.finish_and_clear();

    let count = sequence_lengths.len();

    println!("# records: {count}");
    if count > 0 {
        print_sequence_statistics(
            &mut sequence_lengths,
            &mut sequence_lengths_without_ns,
            additional_percentiles,
            "",
        );
        print_sequence_statistics(
            &mut sequence_hoco_lengths,
            &mut sequence_hoco_lengths_without_ns,
            additional_percentiles,
            "hoco ",
        );
    }

    Ok(())
}

fn print_sequence_statistics(
    sequence_lengths: &mut [usize],
    sequence_lengths_without_ns: &mut [usize],
    additional_percentiles: &[u8],
    prefix: &str,
) {
    sequence_lengths.sort_unstable_by(|a, b| b.cmp(a));
    sequence_lengths_without_ns.sort_unstable_by(|a, b| b.cmp(a));
    let length = sequence_lengths.iter().sum();
    let length_without_ns = sequence_lengths_without_ns.iter().sum();
    let ns = length - length_without_ns;

    println!("{prefix}# Ns: {ns}");
    print_nx(sequence_lengths, length, additional_percentiles, prefix, "");
    print_nx(
        sequence_lengths_without_ns,
        length_without_ns,
        additional_percentiles,
        prefix,
        " without Ns",
    );
}

fn print_nx(
    sorted_sequence_lengths: &[usize],
    length: usize,
    additional_percentiles: &[u8],
    prefix: &str,
    suffix: &str,
) {
    let n50 = nx(sorted_sequence_lengths, length, |l| l / 2);
    let n75 = nx(sorted_sequence_lengths, length, |l| {
        l.checked_mul(3).unwrap() / 4
    });

    println!("{prefix}total length{suffix}: {length}");
    println!("{prefix}N50{suffix}: {n50}");
    println!("{prefix}N75{suffix}: {n75}");

    for additional_percentile in additional_percentiles.iter().copied() {
        let nx = nx(sorted_sequence_lengths, length, |l| {
            ((l as u128) * u128::from(additional_percentile) / 100) as usize
        });
        println!("{prefix}N{additional_percentile}{suffix}: {nx}");
    }

    println!(
        "{prefix}max len{suffix}: {}",
        sorted_sequence_lengths.first().unwrap()
    );
    println!(
        "{prefix}min len{suffix}: {}",
        sorted_sequence_lengths.last().unwrap()
    );
}

fn nx(lengths: &[usize], sum: usize, percentile: impl FnOnce(usize) -> usize) -> usize {
    debug_assert!(lengths.windows(2).all(|w| w[0] >= w[1]));
    debug_assert_eq!(lengths.iter().sum::<usize>(), sum);

    let required_covered_bases = percentile(sum);
    debug_assert!(required_covered_bases <= sum);

    let mut sum = 0;
    for len in lengths.iter().copied() {
        sum += len;
        if sum >= required_covered_bases {
            return len;
        }
    }

    unreachable!()
}

struct SequenceStatistics {
    len: usize,
    hoco_len: usize,
    len_without_ns: usize,
    hoco_len_without_ns: usize,
}

impl SequenceStatistics {
    fn new(sequence: &[u8]) -> Self {
        if sequence.is_empty() {
            return Self {
                len: 0,
                hoco_len: 0,
                len_without_ns: 0,
                hoco_len_without_ns: 0,
            };
        }

        let is_n = |b| b == b'n' || b == b'N';
        let mut len = 1;
        let mut hoco_len = 1;
        let mut last_byte = *sequence.first().unwrap();
        let mut ns = if is_n(last_byte) { 1 } else { 0 };
        let mut hoco_ns = ns;

        for byte in sequence.iter().skip(1).copied() {
            if byte == b'\n' {
                continue;
            }

            len += 1;
            if is_n(byte) {
                ns += 1;
            }

            if byte != last_byte {
                last_byte = byte;
                hoco_len += 1;

                if is_n(last_byte) {
                    hoco_ns += 1;
                }
            }
        }

        Self {
            len,
            hoco_len,
            len_without_ns: len - ns,
            hoco_len_without_ns: hoco_len - hoco_ns,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::basic_statistics;

    #[test]
    fn test() {
        let fasta = b">1\nAAAGCGCTNNNNNTTCGAGGA\n>2\nGTGCTAGCGGGCC\nNCCCTTTTTTTTTTTT\n>3\nACGCTTATG\n>4\nGCTAACTGAGNNNNAAATTTCGGG\n>5\nAAAGGGCCTTCC\n";
        basic_statistics(fasta.as_slice(), fasta.len() as u64, &[], &[90]).unwrap();
    }
}
