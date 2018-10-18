/* Mock Realtime data server
 *
 * Serve "realtime" GTFS data by request? i.e., use worker's pid to serve the next file for that worker.
 */

const express = require('express')
const app = express()

// load all of the files in the ZIP archieve ...
const fs = require('fs')
const path = require('path')
var all_files = fs.readdirSync('archive')
const files = all_files.filter(file => file.match('^vehicle'))

var pids = [] // {"pid001": 1}

app.use('/archive', express.static('archive'))

app.use('/:pid/:start/:end/vehicle_positions', (req, res) => {
    const pid = req.params.pid
    const start = parseInt(req.params.start)
    const end = parseInt(req.params.end)
    if (!(pid in pids)) pids[pid] = start
    if (pids[pid] < start) {
        pids[pid] = start
    }
    if (pids[pid] >= end) {
        res.send('simulation complete')
        return
    }
    res.download(path.join(__dirname, 'archive/' + files[pids[pid]]));
    pids[pid]++
})

app.use('/:pid/vehicle_positions', (req, res) => {
    const pid = req.params.pid
    if (!(pid in pids)) pids[pid] = 0
    if (pids[pid] >= files.length) {
        res.send('simulation complete')
        return
    }
    res.download(path.join(__dirname, 'archive/' + files[pids[pid]]));
    pids[pid]++

    // console.log(pids)
})
app.use('/:pid/reset', (req, res) => {
    const pid = req.params.pid
    if (pid in pids) pids[pid] = 0
    res.send('ok')
})

app.get('*', (req, res) => {
    console.log(req)
    res.status(404).send('Not found')
})

app.listen(3000, () => console.log('Mock GTFS server running on port 3000!'))
