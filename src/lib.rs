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
    input.into_iter().scan(None, |previous_item, item| {
        if let Some(previous_item) = previous_item.as_mut() {
            if *previous_item == item {
                None
            } else {
                *previous_item = item.clone();
                Some(item)
            }
        } else {
            *previous_item = Some(item.clone());
            Some(item)
        }
    })
}
