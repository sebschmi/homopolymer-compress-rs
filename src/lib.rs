//! Functions to homopolymer compress arbitrary sequences.

#![warn(missing_docs)]

/// Homopolymer compress the given sequence.
pub fn homopolymer_compress<
    'output,
    Input: 'output + IntoIterator<Item = Item>,
    Item: 'output + Eq + Clone,
>(
    input: Input,
) -> impl 'output + Iterator<Item = Item> {
    input
        .into_iter()
        .scan(None, |previous_item, item| {
            if let Some(previous_item) = previous_item.as_mut() {
                if *previous_item == item {
                    Some(None)
                } else {
                    *previous_item = item.clone();
                    Some(Some(item))
                }
            } else {
                *previous_item = Some(item.clone());
                Some(Some(item))
            }
        })
        .flatten()
}

/// Homopolymer compress the given sequence and compute a map to homopolymer decompress the output.
pub fn homopolymer_compress_with_hodeco_map<
    'output,
    Input: 'output + IntoIterator<Item = Item>,
    Item: 'output + Eq + Clone,
>(
    input: Input,
) -> impl 'output + Iterator<Item = (Item, usize)> {
    input
        .into_iter()
        .enumerate()
        .scan(None, |previous_item, (index, item)| {
            if let Some(previous_item) = previous_item.as_mut() {
                if *previous_item == item {
                    Some(None)
                } else {
                    *previous_item = item.clone();
                    Some(Some((item, index)))
                }
            } else {
                *previous_item = Some(item.clone());
                Some(Some((item, index)))
            }
        })
        .flatten()
}

#[cfg(test)]
mod tests {
    use crate::{homopolymer_compress, homopolymer_compress_with_hodeco_map};
    use std::iter;

    #[test]
    fn test_homopolymer_compression() {
        let input = b"ACAARRRTGGGTGTJASAAAI";
        let expected = Vec::from_iter(b"ACARTGTGTJASAI".iter().cloned());
        let actual = Vec::from_iter(homopolymer_compress(input.iter().cloned()));
        assert_eq!(expected, actual);
    }

    #[test]
    fn test_hodeco_mapping() {
        let input = b"ACAARRRTGGGTGTJASAAAI";
        let expected = Vec::from_iter(b"ACARTGTGTJASAI".iter().cloned());
        let (actual, mut hodeco_mapping): (Vec<_>, Vec<_>) =
            homopolymer_compress_with_hodeco_map(input.iter().cloned()).unzip();
        hodeco_mapping.push(input.len());
        let hodeco_mapping = hodeco_mapping;
        assert_eq!(expected, actual);
        let hodeco: Vec<_> = actual
            .into_iter()
            .zip(hodeco_mapping.windows(2))
            .flat_map(|(item, count)| iter::repeat(item).take(count[1] - count[0]))
            .collect();
        assert_eq!(hodeco, input);
    }
}
