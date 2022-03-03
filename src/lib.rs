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

#[cfg(test)]
mod tests {
    use crate::homopolymer_compress;

    #[test]
    fn test_homopolymer_compression() {
        let input = b"ACAARRRTGGGTGTJASAAAI";
        let expected = Vec::from_iter(b"ACARTGTGTJASAI".iter().cloned());
        let actual = Vec::from_iter(homopolymer_compress(input.iter().cloned()));
        assert_eq!(expected, actual);
    }
}
