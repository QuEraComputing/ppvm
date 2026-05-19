// Section index for the Python Quick Start. Used by every page under
// `quickstart/python/*` so the left-hand section nav stays consistent and
// the order is defined in a single place.
const base = (import.meta.env.BASE_URL || "/").replace(/\/$/, "");

export interface PyQuickStartSection {
  id: string;
  label: string;
  href: string;
}

export const pyQuickStartSections: PyQuickStartSection[] = [
  { id: "install",          label: "§ 1 Install",                    href: `${base}/quickstart/python` },
  { id: "stim",             label: "§ 2 Stim circuits & sampling",   href: `${base}/quickstart/python/stim` },
  { id: "tableau",          label: "§ 3 Generalized Tableau",        href: `${base}/quickstart/python/tableau` },
  { id: "pauli-propagation", label: "§ 4 Pauli Propagation",         href: `${base}/quickstart/python/pauli-propagation` },
  { id: "loss-channel",     label: "§ 5 Loss channel details",       href: `${base}/quickstart/python/loss-channel` },
  { id: "next-steps",       label: "§ 6 Next steps",                 href: `${base}/quickstart/python/next-steps` },
];
