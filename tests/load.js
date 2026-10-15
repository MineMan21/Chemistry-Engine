// Loads the browser scripts into Node (they attach to globalThis.CE).
const path = require('path');
for (const f of ['elements.js', 'chemistry.js', 'appearance.js', 'phase.js', 'sim.js']) require(path.join(__dirname, '..', 'js', f));
module.exports = globalThis.CE;
