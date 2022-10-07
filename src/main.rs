use clap::Parser;
use seq_io::fastx::Reader;
use seq_io::BaseRecord;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;

#[derive(Parser)]
struct Cli {
    /// Fasta or fastq input file (automatically detected).
    #[clap(index = 1)]
    input: PathBuf,
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();

    if !cli.input.is_file() {
        return Err(format!("Not a file: {:?}", cli.input));
    }

    let input_file =
        File::open(&cli.input).map_err(|err| format!("Cannot open input file: {}", err))?;
    basic_statistics(input_file)
}

fn basic_statistics(input: impl Read) -> Result<(), String> {
    let mut fastx_reader = Reader::new(BufReader::new(input));

    let mut sequence_lengths = Vec::new();
    let mut sequence_hoco_lengths = Vec::new();

    while let Some(record) = fastx_reader.next() {
        let record = record.map_err(|err| format!("Error parsing fastx: {}", err))?;

        sequence_lengths.push(record.seq().len());
        sequence_hoco_lengths.push(hoco_len(record.seq()));
    }

    sequence_lengths.sort_unstable_by(|a, b| b.cmp(a));
    sequence_hoco_lengths.sort_unstable_by(|a, b| b.cmp(a));

    let count = sequence_lengths.len();

    let length = sequence_lengths.iter().sum();
    let hoco_length = sequence_hoco_lengths.iter().sum();

    let n50 = nx(&sequence_lengths, length, |l| l / 2);
    let n75 = nx(&sequence_lengths, length, |l| l.checked_mul(3).unwrap() / 4);
    let n50_hoco = nx(&sequence_hoco_lengths, hoco_length, |l| l / 2);
    let n75_hoco = nx(&sequence_hoco_lengths, hoco_length, |l| {
        l.checked_mul(3).unwrap() / 4
    });

    println!("# records: {count}");
    println!("total length: {length}");
    println!("n50: {n50}");
    println!("n75: {n75}");
    println!("total length hoco: {hoco_length}");
    println!("n50 hoco: {n50_hoco}");
    println!("n75 hoco: {n75_hoco}");

    Ok(())
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

fn hoco_len(sequence: &[u8]) -> usize {
    if sequence.is_empty() {
        return 0;
    }

    let mut len = 1;
    let mut last_byte = *sequence.first().unwrap();
    for byte in sequence.iter().skip(1).copied() {
        if byte != last_byte {
            last_byte = byte;
            len += 1;
        }
    }

    len
}

#[cfg(test)]
mod tests {
    use crate::basic_statistics;

    #[test]
    fn test() {
        let fasta = b">1\nAAAGCGCTTTCGAGGA\n>2\nGTGCTAGCGGGCCCCCTTTTTTTTTTTT\n>3\nACGCTTATG\n>4\nGCTAACTGAGAAATTTCGGG\n>5\nAAAGGGCCTTCC\n";
        basic_statistics(fasta.as_slice()).unwrap();
    }
}
